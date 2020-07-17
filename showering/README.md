# Showering


### 1. For showering with PYTHIA8, user should install the following packages at first:

----
PYTHIA8: http://home.thep.lu.se/Pythia/ 

FASTJET: http://fastjet.fr/quickstart.html

FASTJET-contrib: https://fastjet.hepforge.org/contrib/ 

---

### 2. In order to use `Ghost_Associated_B_hadron.C` to shower and preselection

```
source setset.sh
make Ghost_Associated_B_hadron
```
If nothing wrong, then you can do showering and preselection in fthis way

```
./Ghost_Associated_B_hadron.exe process.cmnd ./preselection_dir/preselection_ <path-to-parton-level-file>/event.lhe. ./truth_dir/truth_ > out.log
```
* process.cmnd: ggh_boost.cmnd, vbf_boost.cmnd, vh_boost.cmnd, tth_boost.cmnd
* preselection_dir: the data after preselection will be stored in this directory ```preselection_.csv``` file.
    ```
		// This will result in a .csv file with format:
		// particle pT, eta, phi, m, id, isCharged \n  \
		// particle pT, eta, phi, m, id, isCharged \n  |
		// ...                                         > an event
		// particle pT, eta, phi, m, id, isCharged \n  |
		// higgs jet pT, eta, phi, m.\n                /
		// \n                                            an empty line
		// particle pT, eta, phi, m, id, isCharged \n  \
		// particle pT, eta, phi, m, id, isCharged \n  |
		// ...                                         > an event
		// particle pT, eta, phi, m, id, isCharged \n  |
		// higgs jet pT, eta, phi, m, evtweight \n     /
		// \n                                            an empty line
    ```
* truth_dir: the truth record for two b quarks will be stored in this directory `truth_.csv` file.
     ```
    // This will result in a .csv file with format:
    // b_E, b_Px, b_Py, b_Pz, bbar_E, bbar_Px, bbar_Py, bbar_Pz, evtweight \n > an event
    ```
* out.log: The presesction and Xection after matching will be in the bottom of this file.
           You also can read the particles infomation step by step in this file.
       
### 3. `loop_out.py` can help to read the preselection rate and matching efficiency in your `.log`.
	You should put all `.log` which belong to a production in a directory.
	 
	```
	python loop_out.py <where-is-your-log>/
	```




