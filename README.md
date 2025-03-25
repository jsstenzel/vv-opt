# vv-opt
Unfortunately requires Python 3.8.8 specifically -- some of the libraries start to fight each other in later Python versions, last time I checked

To get this working in my preferred environment, Windows bash, follow these steps:

1. py -m venv .venv
2. source .venv/Scripts/Activate (each time)
3. py -m pip install -U pip
4. py -m pip install -r requirements.txt

Installing numpy==1.25.2 is unfortunately necessary to get pytensor working outside of an Anaconda environment, but thats the source of all of those pytensor.configdefults warning messages.