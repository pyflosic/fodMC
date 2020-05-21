export PYTHONPATH=../src/:$PYTHONPATH
export PYTHONPATH=../external/:$PYTHONPATH
# Database file for fodmc 
cp ../external/xx_database_xx .
python3 run.py 
