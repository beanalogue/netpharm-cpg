# netpharm

한약-한약 / 한약-양약 상호작용 예측을 위한 네트워크약리학 분석 파이프라인

## 목적
국제 in-silico 학술지 게재 수준의 분석 자동화

## 파이프라인 구성
1. **데이터 수집** (`src/collect.py`) - TCMSP, HERB, BATMAN-TCM, DrugBank, DGIdb
2. **네트워크 구축** (`src/network.py`) - networkx + STRING PPI + Cytoscape 파일 생성
3. **농축분석** (`src/enrichment.py`) - Enrichr API (GO/KEGG)
4. **분자도킹** (`docking/`) - AutoDock Vina 소규모 배치 (선택)
5. **LLM 해석** - Ollama (deepseek-r1:8b) 결과 해석 및 논문 초안

## 환경 설정
```bash
cd ~/netpharm
source venv/bin/activate
```

## 디렉토리 구조
```
netpharm/
├── venv/           # Python 3.12 가상환경
├── data/
│   ├── raw/        # 원본 DB 캐시
│   └── processed/  # 정제 데이터
├── db/             # SQLite 로컬 캐시
├── networks/       # .graphml, Cytoscape 파일
├── enrichment/     # GO/KEGG 결과
├── docking/        # Vina 입출력
├── results/        # 최종 분석 결과
├── notebooks/      # Jupyter 탐색
└── src/            # 파이프라인 소스
```

## 환경 정보
- Server: Intel N150, 11GiB RAM (venv 기반)
- Python: 3.12.3
- LLM: Ollama deepseek-r1:8b (로컬)
