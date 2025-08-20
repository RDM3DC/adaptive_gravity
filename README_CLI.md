# CLI quick start

Install (editable dev):
```bash
pip install -r requirements.txt
pip install -e .
```

Run a galaxy:
```bash
python -c "from adaptive_gravity.simulate import run; print(run('galaxies/NGC_2403.json','runs/NGC_2403'))"
```

Outputs:
- runs/<name>/profile.csv
- runs/<name>/rotation_curve.png
