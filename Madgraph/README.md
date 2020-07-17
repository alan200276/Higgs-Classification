# MadGraph5

#### version: MadGraph5_aMC@NLO.v2.7.2

---
### ggF
1. Download the `The Higgs Characterisation model` and put it in to your madgraph's models file.
>The Higgs Characterisation model: https://feynrules.irmp.ucl.ac.be/wiki/HiggsCharacterisation

```
wget http://feynrules.irmp.ucl.ac.be/raw-attachment/wiki/HiggsCharacterisation/HC_NLO_X0_UFO.zip
```
 then use the folloing comment to excute and generate ggF samples.
 
 ```
 python <where-is-your-madgraph>/bin/mg5_aMC ggh.txt
 ```
 
 * the details are in `ggh.txt`
 
 
 ### VBF, Vh and ttH
1.Using the folloing comment to excute and generate ggF samples.
 
 ```
 python <where-is-your-madgraph>/bin/mg5_aMC <vbf.txt or vh.txt or tth.txt
 ```
 
 * the details are in `.txt` file
 
 
 
 
 
