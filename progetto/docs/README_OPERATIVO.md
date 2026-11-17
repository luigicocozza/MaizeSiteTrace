# README operativo

## Comandi principali

cd ~/progetti/ncrna_rnacentral_ai
source .venv/bin/activate
./run_project.sh
python 17_validate_end_to_end.py --verbose
cat outputs/e2e/e2e_validation_summary.txt
./stop_project.sh
