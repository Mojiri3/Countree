# Countree

<img src="https://github.com/user-attachments/assets/2f862508-dd50-4342-ab0f-2db6beac77ce" height="400"/>

taxid의 출현 빈도에 따라 tree를 가시적으로 보여줌. html로 하면 특정 노드의 정보를 저장할 수도 있음. 대박이지

# 어케 씀
1. ```countree.py``` 다운로드 받아서 실행하셈

# 입력 parameter 뭐임
- 1, 2번째 parameter로는 각각 input 파일과 output 파일을 기입하면 됨. input은 1번째 field에 taxid, 2번째 field에 count 또는 ratio 정보를 지닌 tsv 파일임. output은 png, svg, pdf, html을 확장자로 지님
- --taxdump로 경로 넣으면(안 넣으면 다운 받아서 씀) local db 기반으로 쓸 수 있음
- --layout을 통해 circular, linear을 설정 가능(default: circular인데 솔직히 linear가 보기 더 좋음 키키)
- --info를 통해 input의 2번째 field에 넣은게 count인지 ratio인지 정할 수 있음 (default: ratio)
