export PYTHONPATH=../fodMC/gui/src/:$PYTHONPATH
export PYTHONPATH=../fodMC/gui/external/:$PYTHONPATH
# Database file for fodmc 
cp ../fodMC/gui/external/xx_database_xx .
python3 run.py 
