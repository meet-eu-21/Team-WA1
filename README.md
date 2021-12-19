# Team-WA1

Meet-EU Team WA1

Topic A : Prediction of TADs

## How to download repo? (Linux / MacOS)
1. Create GitHub account
2. Provide email to get access.
3. Create ssh key
* `cd /home/my-user/.ssh`
* `ssh-keygen -t rsa -b 4096`
* provide name and optional password
* `cat <given-name>.pub`
* copy retrieved
* GitHub settings
* SSH and GCP keys tab
* add ssh key
* type any name
* paste result of `cat`-a
* navigate to terminal
* `ssh -T git@github.com -i ~/.ssh/<given-name>`
* now I can:
* `git clone git@github.com:<scope-user>/<project-name>.git`, in our case `git clone git@github.com:kot-sebastian/ipz-tad.git`

## Data description:
1. HiC - data to our algorithms
2. TAD - additional metadata + results

## Downloading data:
1. All (approx. 25 gb) - `./scripts/download_all.sh`
2. Only GM12878 (approx. 11 gb) - `./scripts/download_GM12878.sh`

## Run:
1. install `python3`
2. install required packages
3. download and unpack data or adjust paths in script (file `main.py`)
4. in src directory `python3 main.py`

## Results:
1. Results should be available in `results/topdom` directory (arrowhead is disabled)
2. Currently, script is tailored to 25kb data and results in directory base on this data (you can use 100kb, but result names would not adjust)

