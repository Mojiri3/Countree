import sys
import argparse
import os
import tarfile
import sqlite3
from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, TextFace
import matplotlib
matplotlib.use('Agg')
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import json
import traceback

class CustomNCBITaxa:
    def __init__(self, db_path):
        self.db_path = db_path

    def get_connection(self):
        return sqlite3.connect(self.db_path)

    def get_taxid_translator(self, taxids):
        with self.get_connection() as conn:
            cursor = conn.cursor()
            placeholders = ','.join('?' * len(taxids))
            query = f"SELECT taxid, name FROM names WHERE taxid IN ({placeholders}) AND name_class = 'scientific name'"
            cursor.execute(query, taxids)
            return dict(cursor.fetchall())

    def get_lineage(self, taxid):
        with self.get_connection() as conn:
            cursor = conn.cursor()
            lineage = []
            while taxid != 1:
                cursor.execute("SELECT parent FROM nodes WHERE taxid = ?", (taxid,))
                parent = cursor.fetchone()
                if parent is None:
                    break
                lineage.append(taxid)
                taxid = parent[0]
            return lineage[::-1]

    def get_rank(self, taxid):
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT rank FROM nodes WHERE taxid = ?", (taxid,))
            result = cursor.fetchone()
            return result[0] if result else "Unknown"

    def get_topology(self, taxids):
        common_lineage = None
        for taxid in taxids:
            lineage = self.get_lineage(taxid)
            if common_lineage is None:
                common_lineage = set(lineage)
            else:
                common_lineage.intersection_update(lineage)
        
        tree = Tree()
        for taxid in taxids:
            lineage = self.get_lineage(taxid)
            current_node = tree
            for ancestor in lineage:
                if ancestor not in common_lineage:
                    child = next((child for child in current_node.children if child.name == str(ancestor)), None)
                    if child is None:
                        rank = self.get_rank(ancestor)
                        current_node = current_node.add_child(name=str(ancestor))
                        current_node.add_feature("rank", rank)
                        current_node.add_feature("sci_name", self.get_taxid_translator([ancestor])[ancestor])
                    else:
                        current_node = child
        return tree

def process_taxdump(taxdump_path):
    temp_dir = "temp_taxdump"
    os.makedirs(temp_dir, exist_ok=True)
    
    db_path = os.path.join(temp_dir, "taxa.sqlite")
    
    if os.path.exists(db_path):
        conn = sqlite3.connect(db_path)
        c = conn.cursor()
        c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='nodes'")
        if c.fetchone():
            c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='names'")
            if c.fetchone():
                print("Using existing database.")
                conn.close()
                return db_path

    with tarfile.open(taxdump_path, "r:gz") as tar:
        tar.extractall(path=temp_dir)
    
    nodes_path = os.path.join(temp_dir, "nodes.dmp")
    names_path = os.path.join(temp_dir, "names.dmp")
    
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    
    c.execute('''CREATE TABLE IF NOT EXISTS nodes
                 (taxid INTEGER PRIMARY KEY, parent INTEGER, rank TEXT)''')
    c.execute('''CREATE TABLE IF NOT EXISTS names
                 (taxid INTEGER, name TEXT, name_class TEXT)''')
    
    c.execute("DELETE FROM nodes")
    c.execute("DELETE FROM names")
    
    conn.commit()
    
    with open(nodes_path, 'r') as f:
        for line in f:
            fields = line.strip().split('|')
            taxid = int(fields[0])
            parent = int(fields[1])
            rank = fields[2].strip()
            c.execute("INSERT INTO nodes VALUES (?, ?, ?)", (taxid, parent, rank))
    
    with open(names_path, 'r') as f:
        for line in f:
            fields = line.strip().split('|')
            taxid = int(fields[0])
            name = fields[1].strip()
            name_class = fields[3].strip()
            c.execute("INSERT INTO names VALUES (?, ?, ?)", (taxid, name, name_class))
    
    conn.commit()
    conn.close()
    
    print("Database created or updated successfully.")
    return db_path

def read_data(filename, info_type):
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    taxid = int(parts[0])
                    count = float(parts[1])
                    whole_count = float(parts[2])
                    if info_type == 'count':
                        value = count
                        data[taxid] = {'value': value, 'count': count}
                    else:  # ratio
                        value = count / whole_count if whole_count != 0 else 0
                        data[taxid] = {'value': value, 'count': count, 'whole_count': whole_count}
                except ValueError:
                    print(f"Warning: Skipping invalid line: {line.strip()}")
    return data

def get_threshold_indices(data_values):
    sorted_values = sorted(data_values, reverse=True)
    total = sum(sorted_values)
    cumulative = 0
    thresholds = []
    for i, value in enumerate(sorted_values):
        cumulative += value
        if cumulative >= total * (1 - 1/2**len(thresholds)):
            thresholds.append(value)
        if len(thresholds) == 7:
            break
    return thresholds

def accumulate_values(node, info_type):
    if not node.children:
        if info_type == 'count':
            return node.count, node.value
        else:
            return node.count, node.whole_count, node.value

    total_count = node.count
    total_value = node.value
    if info_type == 'ratio':
        total_whole_count = node.whole_count

    for child in node.children:
        if info_type == 'count':
            child_count, child_value = accumulate_values(child, info_type)
            total_count += child_count
            total_value += child_value
        else:
            child_count, child_whole_count, child_value = accumulate_values(child, info_type)
            total_count += child_count
            total_whole_count += child_whole_count

    if info_type == 'ratio':
        total_value = total_count / total_whole_count if total_whole_count != 0 else 0

    node.add_feature("cumulative_count", total_count)
    node.add_feature("cumulative_value", total_value)
    if info_type == 'ratio':
        node.add_feature("cumulative_whole_count", total_whole_count)
        return total_count, total_whole_count, total_value
    else:
        return total_count, total_value

def create_tree(data, ncbi, layout, info_type):
    valid_taxids = [taxid for taxid in data.keys() if ncbi.get_taxid_translator([taxid])]

    if not valid_taxids:
        raise ValueError("No valid taxids found in the input data.")

    tree = ncbi.get_topology(valid_taxids)

    ts = TreeStyle()
    ts.mode = "c" if layout == "circular" else "r"
    ts.show_leaf_name = True
    ts.scale = 20
    ts.branch_vertical_margin = 10
    ts.show_branch_support = False

    max_count = max(d['count'] for d in data.values())
    threshold_indices = get_threshold_indices([d['value'] for d in data.values()])

    for node in tree.traverse():
        process_node(node, data, ncbi, max_count, threshold_indices, info_type)

    accumulate_values(tree, info_type)

    return tree, ts

def process_node(node, data, ncbi, max_count, threshold_indices, info_type):
    if not node.name:
        node.name = "Unknown"
        taxid = 0
    else:
        try:
            taxid = int(node.name)
        except ValueError:
            print(f"Warning: Invalid taxid '{node.name}'. Setting to 0.")
            taxid = 0
    
    node_data = data.get(taxid, {'count': 0, 'value': 0})
    if info_type == 'ratio':
        node_data.setdefault('whole_count', 0)
    
    count = node_data['count']
    value = node_data['value']
    
    nstyle = NodeStyle()
    nstyle["size"] = 0
    
    style_index = next((i for i, threshold in enumerate(threshold_indices) if value >= threshold), len(threshold_indices))
    
    width = int((style_index + 1) * 2)
    nstyle["hz_line_width"] = max(1, width)
    nstyle["vt_line_width"] = max(1, width)
    
    sci_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
    node.name = f"{sci_name} (ID: {taxid}, Count: {count:.2f}, Value: {value:.4f})"
    node.set_style(nstyle)

    face = CircleFace(radius=max(1, width), color="skyblue", style="sphere")
    node.add_face(face, column=0, position="branch-right")
    
    if hasattr(node, "rank") and hasattr(node, "sci_name"):
        info_text = f"{node.rank}: {node.sci_name}"
        info_face = TextFace(info_text, fsize=8, fgcolor="gray")
        node.add_face(info_face, column=0, position="branch-top")

    node.add_feature("count", count)
    node.add_feature("value", value)
    if info_type == 'ratio':
        node.add_feature("whole_count", node_data['whole_count'])

def create_html_output(tree, output_file, layout, info_type, font, level):
    try:
        print("Starting HTML output creation...")
        tree_data = tree_to_dict(tree)
        print(f"Tree data created. Root node: {tree_data['name']}")
        
        html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Phylogenetic Tree</title>
            <script src="https://d3js.org/d3.v6.min.js"></script>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 0; padding: 0; }}
                #tree-container {{ width: 100vw; height: 100vh; }}
                .node circle {{ fill: #fff; stroke: steelblue; stroke-width: 1.5px; }}
                .node text {{ font-weight: bold; }}
                .link {{ fill: none; stroke: #ccc; stroke-width: 1.5px; }}
                .tooltip {{
                    position: fixed;
                    background-color: white;
                    padding: 10px;
                    border: 1px solid black;
                    border-radius: 5px;
                    opacity: 0;
                    transition: opacity 200ms;
                    pointer-events: auto;
                    user-select: text;
                    z-index: 1000;
                }}
                .context-menu {{
                    position: fixed;
                    background-color: #f9f9f9;
                    padding: 5px;
                    border: 1px solid #ddd;
                    border-radius: 5px;
                    opacity: 0;
                    transition: opacity 200ms;
                    pointer-events: auto;
                    z-index: 1001;
                }}
                .context-menu button {{
                    display: block;
                    width: 100%;
                    padding: 5px;
                    margin: 2px 0;
                    text-align: left;
                    border: none;
                    background: none;
                    cursor: pointer;
                }}
                .context-menu button:hover {{
                    background-color: #f1f1f1;
                }}
                .button-container {{
                    position: fixed;
                    top: 10px;
                    left: 10px;
                    display: flex;
                    flex-direction: column;
                    gap: 10px;
                }}
                .button-container button {{
                    width: 150px;
                }}
                #custom-alert {{
                    position: fixed;
                    color: white;
                    padding: 10px;
                    border-radius: 5px;
                    display: none;
                    z-index: 1000;
                    pointer-events: none;
                }}
                #custom-alert.blue {{
                    background-color: rgba(0, 0, 255, 0.5);
                }}
                #custom-alert.yellow {{
                    background-color: rgba(255, 255, 0, 0.7);
                    color: black;
                }}
                #stored-data-display {{
                    position: fixed;
                    top: 50%;
                    left: 50%;
                    transform: translate(-50%, -50%);
                    background-color: white;
                    padding: 20px;
                    border: 1px solid black;
                    z-index: 1000;
                    max-width: 80%;
                    max-height: 80%;
                    overflow: auto;
                    white-space: pre-wrap;
                    user-select: text;
                    display: none;
                }}
                #legend-container {{
                    position: fixed;
                    top: 20px;
                    right: 280px;
                    width: 200px;
                    font-size: 12px;
                    background-color: rgba(255, 255, 255, 0.8);
                    padding: 10px;
                    border-radius: 5px;
                    box-shadow: 0 0 10px rgba(0,0,0,0.1);
                }}
                #legend-svg {{
                    width: 100%;
                    height: 40px;
                }}
                #node-info-settings {{
                    position: fixed;
                    top: 20px;
                    right: 20px;
                    background-color: white;
                    padding: 10px;
                    border: 1px solid black;
                    z-index: 1000;
                }}
                .node-info-header {{
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                    margin-bottom: 10px;
                }}
                .toggle-button {{
                    background: none;
                    border: none;
                    font-size: 20px;
                    cursor: pointer;
                }}
                .node-info-item {{
                    margin-bottom: 5px;
                }}
                #node-layout-input {{
                    width: 100%;
                    margin-top: 10px;
                }}
                #font-size-slider {{
                    width: 100%;
                    margin-top: 10px;
                }}
                #max-font-size-display {{
                    margin-top: 5px;
                }}
            </style>
        </head>
        <body>
            <div id="tree-container"></div>
            <div id="tooltip" class="tooltip" style="opacity: 0;"></div>
            <div id="context-menu" class="context-menu" style="display: none;">
                <button id="add-button">Add</button>
                <button id="remove-button">Remove</button>
            </div>
            <div class="button-container">
                <button onclick="viewStoredData()">View Stored Data</button>
                <button onclick="removeStoredData()">Remove Stored Data</button>
                <button onclick="storeTheData()">Store The Data</button>
                <button onclick="loadTSVFile()">Load TSV File</button>
            </div>
            <div id="legend-container"></div>
            <div id="custom-alert"></div>
            <div id="stored-data-display"></div>
            <div id="node-info-settings">
                <div class="node-info-header">
                    <h3 style="margin: 0;">Node Display Settings</h3>
                    <button id="toggle-widget" class="toggle-button">−</button>
                </div>
                <div id="widget-content">
                    <div id="node-info-checkboxes"></div>
                    <input type="text" id="node-layout-input" placeholder="Enter node layout">
                    <button onclick="updateNodeDisplay()">Update Display</button>
                    <div>
                        <label for="font-size-slider">Max Font Size:</label>
                        <input type="range" id="font-size-slider" min="10" max="80" value="20">
                        <span id="max-font-size-display">20</span>
                    </div>
                </div>
            </div>
            <script>
            const treeData = {json.dumps(tree_data)};
            const layout = "{layout}";
            const info_type = "{info_type}";
            const font = "{font}";
            const level = "{level}";

            const width = window.innerWidth;
            const height = window.innerHeight;
            const margin = {{top: 50, right: 50, bottom: 50, left: 50}};

            let tree, root;

            if (layout === "circular") {{
                const radius = Math.min(width, height) / 2 - Math.max(margin.top, margin.right, margin.bottom, margin.left);
                tree = d3.cluster().size([2 * Math.PI, radius]);
                root = tree(d3.hierarchy(treeData).sort((a, b) => d3.ascending(a.data.name, b.data.name)));
            }} else {{
                const nodeCount = d3.hierarchy(treeData).descendants().length;
                const treeHeight = Math.max(height - margin.top - margin.bottom, nodeCount * 15);
                tree = d3.tree().size([treeHeight, width - margin.left - margin.right]);
                root = tree(d3.hierarchy(treeData));
            }}

            const svg = d3.select("#tree-container").append("svg")
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("viewBox", layout === "circular" ? [-width / 2, -height / 2, width, height] : [0, 0, width, height])
                .style("font", "10px sans-serif")
                .style("user-select", "none");

            const g = svg.append("g");

            if (layout === "linear") {{
                g.attr("transform", `translate(${{margin.left}},${{margin.top}})`);
            }}

            const cumulative_values = root.descendants().map(d => d.data.cumulative_value);
            const minValue = d3.min(cumulative_values);
            const maxValue = d3.max(cumulative_values);

            let colorScale, fontSizeScale;

            function updateScales(maxFontSize) {{
                const minFontSize = 5;  // 최소 폰트 크기
                const fontSizes = [];
                const step = (maxFontSize - minFontSize) / 8;
                
                for (let i = 1; i < 8; i++) {{  // 가장 작은 값을 무시하고 7개의 값만 사용
                    fontSizes.push(minFontSize + step * i);
                }}

                if (level === "continuous") {{
                    if (font === "absolutely") {{
                        colorScale = d3.scaleSequential(d3.interpolateRainbow)
                            .domain([1, 0]);

                        fontSizeScale = d3.scaleQuantile()
                            .domain([0.00001, 1])
                            .range(fontSizes);
                    }} else {{ // relatively
                        colorScale = d3.scaleSequential(d3.interpolateRainbow)
                            .domain([maxValue, minValue]);

                        fontSizeScale = d3.scaleQuantile()
                            .domain([minValue, maxValue])
                            .range(fontSizes);
                    }}
                }} else {{ // discontinuous
                    const colors = ["#8F00FF", "#4B0082", "#0000FF", "#00FF00", "#FFFF00", "#FF7F00", "#FF0000"];
                    if (font === "absolutely") {{
                        const thresholds = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1]; //Arabidopsis의 chloroplast/all protein 비율이 0.0125정도임
                        colorScale = d3.scaleThreshold()
                            .domain(thresholds)
                            .range(colors);

                        fontSizeScale = d3.scaleThreshold()
                            .domain(thresholds)
                            .range(fontSizes);
                    }} else {{ // relatively
                        const step = (maxValue - minValue) / (colors.length - 1);
                        const thresholds = d3.range(colors.length - 1).map(i => minValue + step * (i + 1));
                        
                        colorScale = d3.scaleThreshold()
                            .domain(thresholds)
                            .range(colors);

                        fontSizeScale = d3.scaleQuantile()
                            .domain([minValue, maxValue])
                            .range(fontSizes);
                    }}
                }}
            }}

            // Initialize scales with default max font size
            updateScales(20);

            function getColor(cumulative_value) {{
                return colorScale(cumulative_value);
            }}

            function getFontSize(cumulative_value) {{
                return fontSizeScale(Math.max(0.00001, cumulative_value));
            }}

            const link = g.append("g")
                .attr("fill", "none")
                .attr("stroke", "#555")
                .attr("stroke-opacity", 0.4)
                .attr("stroke-width", 1.5)
                .selectAll("path")
                .data(root.links())
                .join("path")
                .attr("d", layout === "circular" 
                    ? d3.linkRadial().angle(d => d.x).radius(d => d.y)
                    : d3.linkHorizontal().x(d => d.y).y(d => d.x));


            const toggleButton = document.getElementById('toggle-widget');
            const widgetContent = document.getElementById('widget-content');
            let isWidgetOpen = true;

            toggleButton.addEventListener('click', () => {{
                isWidgetOpen = !isWidgetOpen;
                widgetContent.style.display = isWidgetOpen ? 'block' : 'none';
                toggleButton.textContent = isWidgetOpen ? '−' : '+';
            }});
                
            let nodeInfoDisplay = {{
                'sci_name': true,
                'rank': true,
                'taxid': false,
                'count': false,
                'value': false,
                'cumulative_count': false,
                'cumulative_value': false
            }};

            let nodeLayout = "$(sci_name) ($(rank))";

            function createNodeInfoSettings() {{
                const container = document.getElementById('node-info-checkboxes');
                for (const [key, value] of Object.entries(nodeInfoDisplay)) {{
                    const div = document.createElement('div');
                    div.className = 'node-info-item';
                    const checkbox = document.createElement('input');
                    checkbox.type = 'checkbox';
                    checkbox.id = key;
                    checkbox.checked = value;
                    checkbox.onchange = updateCheckboxState;
                    const label = document.createElement('label');
                    label.htmlFor = key;
                    label.textContent = key;
                    div.appendChild(checkbox);
                    div.appendChild(label);
                    container.appendChild(div);
                }}
                document.getElementById('node-layout-input').value = nodeLayout;
            }}

            function updateCheckboxState(event) {{
                const key = event.target.id;
                nodeInfoDisplay[key] = event.target.checked;
                updateNodeLayout();
            }}

            function updateNodeLayout() {{
                let newLayout = "";
                for (const [key, value] of Object.entries(nodeInfoDisplay)) {{
                    if (value) {{
                        newLayout += `$(${{key}}) `;
                    }}
                }}
                nodeLayout = newLayout.trim();
                document.getElementById('node-layout-input').value = nodeLayout;
                updateNodeDisplay();
            }}

            function updateNodeDisplay() {{
                nodeLayout = document.getElementById('node-layout-input').value;
                updateNodeText();
            }}

            function updateNodeText() {{
                node.selectAll('text')
                    .text(d => getNodeLabel(d.data))
                    .attr("fill", d => getColor(d.data.cumulative_value))
                    .style("font-size", d => `${{getFontSize(d.data.cumulative_value)}}px`);
            }}

            function getNodeLabel(data) {{
                let label = nodeLayout;
                const placeholders = {{
                    '$(sci_name)': data.sci_name,
                    '$(rank)': data.rank,
                    '$(taxid)': data.name.split('(ID:')[1].split(',')[0].trim(),
                    '$(count)': data.count,
                    '$(value)': data.value.toFixed(4),
                    '$(cumulative_count)': data.cumulative_count,
                    '$(cumulative_value)': data.cumulative_value.toFixed(4)
                }};

                for (const [placeholder, value] of Object.entries(placeholders)) {{
                    label = label.replace(placeholder, value);
                }}

                return label.trim();
            }}

            createNodeInfoSettings();

            const node = g.append("g")
                .attr("stroke-linejoin", "round")
                .attr("stroke-width", 3)
                .selectAll("g")
                .data(root.descendants())
                .join("g")
                .attr("transform", d => layout === "circular"
                    ? `rotate(${{d.x * 180 / Math.PI - 90}}) translate(${{d.y}},0)`
                    : `translate(${{d.y}},${{d.x}})`);

            node.append("circle")
                .attr("fill", d => d.children ? "#555" : "#999")
                .attr("r", 2.5);

            node.append("text")
                .attr("dy", "0.31em")
                .attr("x", d => layout === "circular"
                    ? (d.x < Math.PI === !d.children ? 6 : -6)
                    : (d.children ? -6 : 6))
                .attr("text-anchor", d => layout === "circular"
                    ? (d.x < Math.PI === !d.children ? "start" : "end")
                    : (d.children ? "end" : "start"))
                .attr("transform", d => layout === "circular" && d.x >= Math.PI ? "rotate(180)" : null)
                .text(d => getNodeLabel(d.data))
                .attr("fill", d => getColor(d.data.cumulative_value))
                .style("font-size", d => `${{getFontSize(d.data.cumulative_value)}}px`)
                .clone(true).lower()
                .attr("stroke", "white");

            const tooltip = d3.select("#tooltip");
            const contextMenu = d3.select("#context-menu");

            let activeNode = null;
            let hideTimer = null;
            let isMouseOverNode = false;
            let isMouseOverTooltipOrMenu = false;

            function showTooltipAndMenu(event, d) {{
                if (activeNode !== d) {{
                    clearTimeout(hideTimer);
                    hideTooltipAndMenu(() => {{
                        displayNodeInfo(event, d);
                    }});
                }}
                isMouseOverNode = true;
            }}

            function displayNodeInfo(event, d) {{
                const [x, y] = d3.pointer(event, svg.node());
                const adjustedX = x + (layout === "circular" ? width / 2 : 0);
                const adjustedY = y + (layout === "circular" ? height / 2 : 0);
                
                tooltip.style("opacity", 0)
                    .style("left", adjustedX + "px")
                    .style("top", adjustedY + "px")
                    .html(`<div>
                            <strong>${{d.data.sci_name}}</strong><br/>
                            Rank: ${{d.data.rank}}<br/>
                            Taxid: ${{d.data.name.split('(ID:')[1].split(',')[0].trim()}}<br/>
                            Count: ${{d.data.count}}<br/>
                            Value: ${{d.data.value.toFixed(4)}}<br/>
                            Cumulative Count: ${{d.data.cumulative_count}}<br/>
                            Cumulative Value: ${{d.data.cumulative_value.toFixed(4)}}
                            ${{info_type === 'ratio' ? `<br/>Whole Count: ${{d.data.whole_count}}<br/>Cumulative Whole Count: ${{d.data.cumulative_whole_count}}` : ''}}
                           </div>`);

                tooltip.transition().duration(200).style("opacity", 0.9);

                const tooltipRect = tooltip.node().getBoundingClientRect();

                contextMenu.style("opacity", 0)
                    .style("display", "block")
                    .style("left", (tooltipRect.right + 5) + "px")
                    .style("top", tooltipRect.top + "px");

                contextMenu.transition().duration(200).style("opacity", 1);

                d3.select("#add-button").on("click", () => addToTSV(d));
                d3.select("#remove-button").on("click", () => removeFromTSV(d));
                
                activeNode = d;
            }}

            function startHideTimer() {{
                isMouseOverNode = false;
                if (!isMouseOverTooltipOrMenu) {{
                    clearTimeout(hideTimer);
                    hideTimer = setTimeout(() => hideTooltipAndMenu(), 500);
                }}
            }}

            function hideTooltipAndMenu(callback) {{
                if (!activeNode || (!isMouseOverNode && !isMouseOverTooltipOrMenu)) {{
                    tooltip.transition().duration(200).style("opacity", 0);
                    contextMenu.transition().duration(200).style("opacity", 0);
                    setTimeout(() => {{
                        contextMenu.style("display", "none");
                        activeNode = null;
                        if (callback) callback();
                    }}, 200);
                }}
            }}

            node.on("mouseover", (event, d) => showTooltipAndMenu(event, d))
                .on("mouseout", startHideTimer);

            tooltip.on("mouseover", () => {{
                clearTimeout(hideTimer);
                isMouseOverTooltipOrMenu = true;
            }})
            .on("mouseout", () => {{
                isMouseOverTooltipOrMenu = false;
                if (!isMouseOverNode) {{
                    startHideTimer();
                }}
            }});

            contextMenu.on("mouseover", () => {{
                clearTimeout(hideTimer);
                isMouseOverTooltipOrMenu = true;
            }})
            .on("mouseout", () => {{
                isMouseOverTooltipOrMenu = false;
                if (!isMouseOverNode) {{
                    startHideTimer();
                }}
            }});

            function showCustomAlert(message, x, y, isWarning = false) {{
                const alert = document.getElementById('custom-alert');
                alert.style.left = (x + 10) + 'px';
                alert.style.top = (y - 50) + 'px';
                alert.innerHTML = message;
                alert.className = isWarning ? 'yellow' : 'blue';
                alert.style.display = 'block';
                setTimeout(() => {{
                    alert.style.display = 'none';
                }}, 2000);
            }}

            function addToTSV(d) {{
                const taxid = d.data.name.split('(ID:')[1].split(',')[0].trim();
                let stored = JSON.parse(localStorage.getItem('storedData') || '[]');
                const existingIndex = stored.findIndex(item => item.taxid === taxid);
                if (existingIndex === -1) {{
                    stored.push({{
                        rank: d.data.rank,
                        sci_name: d.data.sci_name,
                        taxid: taxid,
                        count: d.data.count,
                        ratio: d.data.ratio,
                        cumulative_count: d.data.cumulative_count,
                        cumulative_value: d.data.cumulative_value
                    }});
                    localStorage.setItem('storedData', JSON.stringify(stored));
                    showCustomAlert(`Added sci_name ${{d.data.sci_name}} (Taxid: ${{taxid}})`, event.pageX, event.pageY);
                }} else {{
                    showCustomAlert(`Taxid ${{taxid}} already exists in stored data`, event.pageX, event.pageY, true);
                }}
            }}

            function removeFromTSV(d) {{
                const taxid = d.data.name.split('(ID:')[1].split(',')[0].trim();
                let stored = JSON.parse(localStorage.getItem('storedData') || '[]');
                stored = stored.filter(item => item.taxid !== taxid);
                localStorage.setItem('storedData', JSON.stringify(stored));
                showCustomAlert(`Removed sci_name ${{d.data.sci_name}} (Taxid: ${{taxid}})`, event.pageX, event.pageY);
            }}

            function viewStoredData() {{
                let stored = JSON.parse(localStorage.getItem('storedData') || '[]');
                const storedDataDisplay = document.getElementById('stored-data-display');
                storedDataDisplay.textContent = JSON.stringify(stored, null, 2);
                storedDataDisplay.style.display = 'block';
                const closeButton = document.createElement('button');
                closeButton.textContent = 'Close';
                closeButton.style.marginTop = '10px';
                closeButton.onclick = () => storedDataDisplay.style.display = 'none';
                
                storedDataDisplay.appendChild(closeButton);
            }}

            function storeTheData() {{
                let stored = JSON.parse(localStorage.getItem('storedData') || '[]');
                if (stored.length === 0) {{
                    alert('No data to store. Please add some data first.');
                    return;
                }}

                const tsvContent = stored.map(item => 
                    `${{item.rank}}\t${{item.sci_name}}\t${{item.taxid}}\t${{item.count}}\t${{item.ratio}}\t${{item.cumulative_count}}\t${{item.cumulative_value}}`
                ).join('\\n');

                const blob = new Blob([tsvContent], {{ type: 'text/tab-separated-values' }});
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                link.href = url;
                const outputFileName = '{os.path.splitext(os.path.basename(output_file))[0]}';
                link.download = `${{outputFileName}}.select`;
                link.click();
                URL.revokeObjectURL(url);
                alert('Data has been stored in the TSV file.');
            }}

            function loadTSVFile() {{
                if (confirm('Are you sure you want to load the TSV file? This will overwrite the current stored data.')) {{
                    const input = document.createElement('input');
                    input.type = 'file';
                    input.accept = '.select';
                    input.onchange = function(event) {{
                        const file = event.target.files[0];
                        const reader = new FileReader();
                        reader.onload = function(e) {{
                            const content = e.target.result;
                            const lines = content.split('\\n');
                            const data = lines.map(line => {{
                                const [rank, sci_name, taxid, count, ratio, cumulative_count, cumulative_value] = line.split('\\t');
                                return {{ 
                                    rank, 
                                    sci_name, 
                                    taxid, 
                                    count: parseFloat(count), 
                                    ratio: parseFloat(ratio),
                                    cumulative_count: parseFloat(cumulative_count),
                                    cumulative_value: parseFloat(cumulative_value)
                                }};
                            }});
                            localStorage.setItem('storedData', JSON.stringify(data));
                            alert('TSV file loaded successfully.');
                        }};
                        reader.readAsText(file);
                    }};
                    input.click();
                }}
            }}

            function removeStoredData() {{
                if (confirm('Are you sure you want to remove all stored data?')) {{
                    localStorage.removeItem('storedData');
                    alert('All stored data has been removed.');
                }}
            }}

            const zoom = d3.zoom()
                .scaleExtent([0.1, 10])
                .on("zoom", (event) => {{
                    g.attr("transform", event.transform);
                    hideTooltipAndMenu();
                    
                    node.selectAll("text")
                        .style("font-size", d => `${{getFontSize(d.data.cumulative_value) / Math.sqrt(event.transform.k)}}px`);
                }});

            svg.call(zoom);

            // Font size slider
            const fontSizeSlider = document.getElementById('font-size-slider');
            const maxFontSizeDisplay = document.getElementById('max-font-size-display');

            fontSizeSlider.addEventListener('input', function() {{
                const maxFontSize = parseInt(this.value);
                maxFontSizeDisplay.textContent = maxFontSize;
                updateScales(maxFontSize);
                updateNodeText();
            }});

            </script>
        </body>
        </html>
        """
        
        print(f"HTML content created. Length: {len(html_content)} characters")

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"HTML file written to {output_file}")
        
        if os.path.exists(output_file):
            print(f"Confirmed: {output_file} exists.")
            print(f"File size: {os.path.getsize(output_file)} bytes")
        else:
            print(f"Error: {output_file} was not created.")
        
    except Exception as e:
        print(f"Error in create_html_output: {str(e)}")
        print("Traceback:")
        traceback.print_exc()

def tree_to_dict(node):
    node_dict = {
        "name": node.name,
        "children": [tree_to_dict(child) for child in node.children],
        "rank": getattr(node, "rank", ""),
        "sci_name": getattr(node, "sci_name", ""),
        "count": getattr(node, "count", 0),
        "value": getattr(node, "value", 0), 
        "cumulative_count": getattr(node, "cumulative_count", getattr(node, "count", 0)),
        "cumulative_value": getattr(node, "cumulative_value", getattr(node, "value", 0))
    }
    if hasattr(node, "whole_count"):
        node_dict["whole_count"] = node.whole_count
    if hasattr(node, "cumulative_whole_count"):
        node_dict["cumulative_whole_count"] = node.cumulative_whole_count
    return node_dict
def main():
    parser = argparse.ArgumentParser(description="Generate a phylogenetic tree from NCBI taxids (TSV input).")
    parser.add_argument("input_file", help="Input TSV file containing taxids and counts/ratios")
    parser.add_argument("output_file", help="Output file name (e.g., tree.html)")
    parser.add_argument("--taxdump", help="Path to taxdump.tar.gz file", required=True)
    parser.add_argument("--layout", choices=["circular", "linear"], default="circular", help="Tree layout (default: circular)")
    parser.add_argument("--info", choices=["count", "ratio"], default="count", help="Type of information in the input file (default: count)")
    parser.add_argument("--font", choices=["absolutely", "relatively"], default="absolutely", help="Font scaling method (default: absolutely)")
    parser.add_argument("--level", choices=["continuous", "discontinuous"], default="discontinuous", help="Color and font scaling method (default: discontinuous)")
    args = parser.parse_args()

    print("Processing taxdump file...")
    db_path = process_taxdump(args.taxdump)
    print("Taxdump processed and database ready.")

    ncbi = CustomNCBITaxa(db_path)

    data = read_data(args.input_file, args.info)
    if not data:
        print("Error: No valid data found in the input file.")
        sys.exit(1)

    if not args.output_file.lower().endswith('.html'):
        print("Error: Output file must be an HTML file.")
        sys.exit(1)

    try:
        print(f"Creating {args.layout} tree...")
        tree, ts = create_tree(data, ncbi, args.layout, args.info)
        print(f"{args.layout.capitalize()} tree created successfully.")

        print("Generating HTML output...")
        create_html_output(tree, args.output_file, args.layout, args.info, args.font, args.level)

        print(f"{args.layout.capitalize()} tree has been saved as '{args.output_file}'")
        print(f"You can open this file in a web browser to view the tree.")

    except Exception as e:
        print(f"Error in main: {str(e)}")
        print("Traceback:")
        traceback.print_exc()

if __name__ == "__main__":
    main()