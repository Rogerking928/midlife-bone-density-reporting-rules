# Postmenopausal and premenopausal bone density reporting rules in midlife women (NHANES)

Analysis code for the cross-classification of femoral neck DXA results in US women aged 45-55 years under the
postmenopausal rule (T-score, WHO categories) and the premenopausal rule (age-, sex- and race-matched Z-score, with -2.0
or lower below the expected range for age), using the National Health and Nutrition Examination Survey (NHANES)
2005-2006, 2007-2008, 2009-2010, 2013-2014 and 2017-2018.

## Running

```
pip install -r requirements.txt
python3 fetch_data.py            # downloads the 35 NHANES files into data/ and checks their MD5
cd code
python3 analysis.py              # base estimates            -> out/results.json
python3 analysis2.py             # sites, splines, fracture  -> out/results2.json
python3 analysis3.py             # counts, T vs Z, variance  -> out/results3.json
python3 analysis4.py             # matched cut-points, sensitivity -> out/results4.json
python3 analysis5.py             # Table 1 extra columns, fracture at 40+ -> out/results5.json
python3 figs.py                  # Figure S1 (participant flow)
python3 figs3.py                 # Figures 1-3 and S2-S5 (figs2.py holds shared drawing helpers)
python3 book3.py                 # Tables_and_Figures.docx, PNGs and Source_data.xlsx -> $DELIVERABLE_DIR
```

`DELIVERABLE_DIR` defaults to `deliverable/` in the repository. All bootstrap steps use fixed seeds. Comments in the
code are partly in Chinese.

## Data

NHANES public-use files are produced by the National Center for Health Statistics and are available at
https://wwwn.cdc.gov/nchs/nhanes/. They are not redistributed here; `fetch_data.py` downloads them.

## Software

Python 3.12.3; package versions in `requirements.txt`.

The analysis was written and run by the author in Python 3, with Claude (Anthropic) and ChatGPT (OpenAI) used as
assistants under the author's direction; the author reviewed and tested all code and takes full responsibility for
the results it produces.

## Licence

MIT.
