# Countree

<img src="https://github.com/user-attachments/assets/372e2fbd-0378-41e0-97eb-7aece971798a" height="400"/>

taxid의 출현 빈도에 따라 tree를 가시적으로 보여줌. html로 하면 특정 노드의 정보를 저장할 수도 있음. 대박이지

# 어케 씀
1. ```countree.py``` 다운로드 받으셈
2. 대부분의 기능은 아래 parameter 설명 잘 읽고 따라해보셈

# 입력 parameter 뭐임
- 1, 2번째 parameter로는 각각 ```input``` 파일과 ```output``` 파일을 기입하면 됨. input은 1번째 field에 taxid, 2번째 field에 count 정보를 지닌 tsv 파일임.
  만약, ```ratio```를 사용할 거라면 3번째 field에 해당 taxid의 전체 개수를 자신의 db에서 긁어오면 됨.
  - ```output```은 png, svg, pdf, html을 확장자로 지님
- ```--taxdump```로 경로 넣으면(안 넣으면 다운 받아서 씀) local db 기반으로 쓸 수 있음
- ```--layout```을 통해 circular, linear을 설정 가능(default: circular인데 솔직히 linear가 보기 더 좋음 키키)
- ```--info```를 통해 input의 2번째 field에 넣은게 count인지 ratio인지 정할 수 있음 (default: ratio)
- ```--font```를 통해 font 서식 차이를 전체 범위(전체 비중의 0\~1 사이 값) 내에서 다룰 것인지, 아니면 최소\~최대 범위에서 다룰 것인지 고름 (default: absolutely)
- ```--level```을 통해 font 서식의 연속성과 불연속성을 정함 (default: discontinuous)

# Input file 어케 만듦
뭐, 방식은 많겠지만, 나는 아래와 같이 함
1. BLAST output을 통해 taxid의 counting 진행
-  나도 뉴비라 보고 쓰려고 커맨드 예시 보여드림

    ```awk -F'\t' -v m=(taxid field 번호) -v n=(subject field 번호) '{a[$m][($n)] = 1} END {for (i in a) print i "\t" length(a[i])}' OFS='\t' (tabular type의 BLAST output) | sort -k2 -nr > (Input file)```
2. ```--info```가 ```count```면 여기서 끝내도 되는데, ```ratio```라면 추가 작업 필요함.
-  3번째 field를 채우기 위해 여러 방법이 있지만, 나는 이렇게 했음.

    ```awk 'BEGIN {FS=OFS="\t"} {count = "esearch -db protein -query '\''txid" $1 "[Organism:noexp] & (\"0001/01/01\"[PDAT] :\"yyyy/mm/dd\"[PDAT]) & refseq[filter]'\'' | xtract -pattern ENTREZ_DIRECT -element Count"; count | getline result; close(count); print $0, result}' (1번에서 제작한 파일) > (Input file)```

# 기대 기능
1. 솔직히 png, svg, pdf는 비추천함. html에 신경 더 써썼으니 그걸로 하셈. 용량도 별로 안 큼.
2. taxid에 해당하는 대상이 얼마나 차지하는지 tree에서 가시적으로 볼 수 있음.
3. output format 이 html일 경우
  - tree에서 node 정보를 즉각적으로 볼 수 있고, data를 저장할 수 있음.
  - 저장한 data를 다운로드 받아 tsv 파일로 받을 수 있음.
  - tsv 파일을 수정한 뒤, upload해서 데이터를 유동적으로 관리 가능함.
4. 생성된 html에 node에 표시되는 text를 바꾸는 기능과 font 크기 조절 기능이 있음. 적용하고 안 된 거 같다 싶으면 확대/축소 ㄱ
## Tree 예시 보여줄게
- ```--info ratio --layout linear --font absolutely --level discontinuous```
<img src="https://github.com/user-attachments/assets/4fb9dbda-55fc-4158-882c-5b4880ccdc61" height="400"/>

- ```--info ratio --layout linear --font relatively --level discontinuous```
<img src="https://github.com/user-attachments/assets/0b86dc4c-a8e6-457a-ae12-ad57e7c48842" height="400"/>

난 ```linear```랑 ```discontinuous``` 를 쓰는 게 예뻐서 ```--font```만 다르게 해서 보여줬음. 딴 건 직접 해봐.

# 한계
1. html로 저장하는 파일이 html 경로에 저장되지 않음. Linux에서 만든 거 ssh 통해서 Windows 환경에서 여는데 계속 download 경로로 받아지더라. 이건 좀 귀찮을 듯.
2. Input file 만드는 게 좀 귀찮긴 함. 근데 어쩔티비????????? 굴릴 때마다 파일 입력하는 거는 작업 시간이 길어지는 문제도 있고, 무엇보다 굳이 불필요한 기능같음. 나중에 입력하기 귀찮으면 만들래.
