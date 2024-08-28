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
    total_count = 0
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    taxid = int(parts[0])
                    value = float(parts[1])
                    data[taxid] = value
                    if info_type == 'count':
                        total_count += value
                except ValueError:
                    print(f"Warning: Skipping invalid line: {line.strip()}")
    
    if info_type == 'count':
        # Calculate ratio for each taxid
        for taxid in data:
            ratio = data[taxid] / total_count
            data[taxid] = {'count': data[taxid], 'ratio': ratio}
    else:
        # For 'ratio' info_type, set count same as value and ratio as value
        for taxid in data:
            data[taxid] = {'count': data[taxid], 'ratio': data[taxid]}
    
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
    threshold_indices = get_threshold_indices([d['ratio'] for d in data.values()])

    for node in tree.traverse():
        process_node(node, data, ncbi, max_count, threshold_indices, info_type)

    accumulate_counts(tree)

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
    
    node_data = data.get(taxid, {'count': 0, 'ratio': 0})
    count = node_data['count']
    ratio = node_data['ratio']
    
    nstyle = NodeStyle()
    nstyle["size"] = 0
    
    # Determine the index for styling based on ratio
    style_index = next((i for i, threshold in enumerate(threshold_indices) if ratio >= threshold), len(threshold_indices))
    
    width = int((style_index + 1) * 2)  # Increase width based on index
    nstyle["hz_line_width"] = max(1, width)
    nstyle["vt_line_width"] = max(1, width)
    
    sci_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
    node.name = f"{sci_name} (ID: {taxid}, Count: {count:.2f}, Ratio: {ratio:.4f})"
    node.set_style(nstyle)

    face = CircleFace(radius=max(1, width), color="skyblue", style="sphere")
    node.add_face(face, column=0, position="branch-right")
    
    if hasattr(node, "rank") and hasattr(node, "sci_name"):
        info_text = f"{node.rank}: {node.sci_name}"
        info_face = TextFace(info_text, fsize=8, fgcolor="gray")
        node.add_face(info_face, column=0, position="branch-top")

    node.add_feature("count", count)
    node.add_feature("ratio", ratio)

def accumulate_counts(node):
    if not hasattr(node, 'count'):
        node.add_feature('count', 0)
    if not hasattr(node, 'ratio'):
        node.add_feature('ratio', 0)

    if not node.children:
        return node.count, node.ratio
    
    total_count = node.count
    total_ratio = node.ratio
    for child in node.children:
        child_count, child_ratio = accumulate_counts(child)
        total_count += child_count
        total_ratio += child_ratio
    
    node.add_feature("cumulative_count", total_count)
    node.add_feature("cumulative_ratio", total_ratio)
    return total_count, total_ratio

def create_html_output(tree, output_file, layout):
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
            <div id="custom-alert"></div>
            <div id="stored-data-display"></div>
            <script>
            const treeData = {json.dumps(tree_data)};
            const layout = "{layout}";

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
                .attr("width", width)
                .attr("height", height)
                .attr("viewBox", layout === "circular" ? [-width / 2, -height / 2, width, height] : [0, 0, width, height])
                .style("font", "10px sans-serif")
                .style("user-select", "none");

            const g = svg.append("g");

            if (layout === "linear") {{
                g.attr("transform", `translate(${{margin.left}},${{margin.top}})`);
            }}
            
            const colorScale = cumulative_ratio => {{
                if (cumulative_ratio >= 0.5) return "#FF0000";  // Red
                if (cumulative_ratio >= 0.25) return "#FF7F00";   // Orange
                if (cumulative_ratio >= 0.125) return "#FFFF00";   // Yellow
                if (cumulative_ratio >= 0.0625) return "#00FF00";   // Green
                if (cumulative_ratio >= 0.03125) return "#0000FF";    // Blue
                if (cumulative_ratio >= 0.015625) return "#4B0082";    // Indigo
                return "#8F00FF";                     // Violet
            }};

            const fontSize = 30;
            
            const fontSizeScale = cumulative_ratio => {{
                if (cumulative_ratio >= 0.5) return fontSize;
                if (cumulative_ratio >= 0.25) return fontSize * 0.9;
                if (cumulative_ratio >= 0.125) return fontSize * 0.8;
                if (cumulative_ratio >= 0.0625) return fontSize * 0.7;
                if (cumulative_ratio >= 0.03125) return fontSize * 0.6;
                if (cumulative_ratio >= 0.015625) return fontSize * 0.5;
                return fontSize * 0.4;
            }};

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
                .text(d => d.data.rank)
                .attr("fill", d => colorScale(d.data.cumulative_ratio))
                .style("font-size", d => `${{fontSizeScale(d.data.cumulative_ratio)}}px`)
                .clone(true).lower()
                .attr("stroke", "white");

            const tooltip = d3.select("#tooltip");
            const contextMenu = d3.select("#context-menu");

            let activeNode = null;
            let hideTimer = null;
            let isHiding = false;

            function showTooltipAndMenu(event, d) {{
                if (isHiding) {{
                    return;
                }}

                if (activeNode !== d) {{
                    clearTimeout(hideTimer);
                    hideTooltipAndMenu(() => {{
                        displayNodeInfo(event, d);
                    }});
                }}
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
                            Ratio: ${{d.data.ratio.toFixed(4)}}<br/>
                            Cumulative Count: ${{d.data.cumulative_count}}<br/>
                            Cumulative Ratio: ${{d.data.cumulative_ratio.toFixed(4)}}
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
                isHiding = false;
            }}

            function startHideTimer() {{
                clearTimeout(hideTimer);
                hideTimer = setTimeout(() => hideTooltipAndMenu(), 500);
            }}

            function hideTooltipAndMenu(callback) {{
                if (!activeNode) {{
                    if (callback) callback();
                    return;
                }}
                
                isHiding = true;

                tooltip.transition().duration(200).style("opacity", 0);
                contextMenu.transition().duration(200).style("opacity", 0);
                setTimeout(() => {{
                    contextMenu.style("display", "none");
                    activeNode = null;
                    isHiding = false;
                    if (callback) callback();
                }}, 200);
            }}

            node.on("mouseover", showTooltipAndMenu)
                .on("mouseout", startHideTimer);

            tooltip.on("mouseover", () => clearTimeout(hideTimer))
                   .on("mouseout", startHideTimer);

            contextMenu.on("mouseover", () => clearTimeout(hideTimer))
                       .on("mouseout", startHideTimer);

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
                        cumulative_ratio: d.data.cumulative_ratio
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
                    `${{item.rank}}\t${{item.sci_name}}\t${{item.taxid}}\t${{item.count}}\t${{item.ratio}}\t${{item.cumulative_count}}\t${{item.cumulative_ratio}}`
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
                                const [rank, sci_name, taxid, count, ratio, cumulative_count, cumulative_ratio] = line.split('\\t');
                                return {{ 
                                    rank, 
                                    sci_name, 
                                    taxid, 
                                    count: parseFloat(count), 
                                    ratio: parseFloat(ratio),
                                    cumulative_count: parseFloat(cumulative_count),
                                    cumulative_ratio: parseFloat(cumulative_ratio)
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
                        .style("font-size", d => `${{fontSizeScale(d.data.cumulative_ratio) / Math.sqrt(event.transform.k)}}px`);
                }});

            svg.call(zoom);

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
        "ratio": getattr(node, "ratio", 0),
        "cumulative_count": getattr(node, "cumulative_count", getattr(node, "count", 0)),
        "cumulative_ratio": getattr(node, "cumulative_ratio", getattr(node, "ratio", 0))
    }
    return node_dict

def main():
    parser = argparse.ArgumentParser(description="Generate a phylogenetic tree from NCBI taxids (TSV input).")
    parser.add_argument("input_file", help="Input TSV file containing taxids and counts/ratios")
    parser.add_argument("output_file", help="Output file name (e.g., tree.html)")
    parser.add_argument("--taxdump", help="Path to taxdump.tar.gz file", required=True)
    parser.add_argument("--layout", choices=["circular", "linear"], default="circular", help="Tree layout (default: circular)")
    parser.add_argument("--info", choices=["count", "ratio"], default="count", help="Type of information in the input file (default: count)")
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
        create_html_output(tree, args.output_file, args.layout)

        print(f"{args.layout.capitalize()} tree has been saved as '{args.output_file}'")
        print(f"You can open this file in a web browser to view the tree.")

    except Exception as e:
        print(f"Error in main: {str(e)}")
        print("Traceback:")
        traceback.print_exc()

if __name__ == "__main__":
    main()