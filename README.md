# cotw

### 🪸🪸 🌏 COTW: Corals Of The World 🌏 🪸🪸
Code and scripts for: A Global Core Coral Metabolome as a Baseline for Holobiont Physiology and Stress Biology  

## Repository structure

```text
cotw/
├── analysis/
│   ├── collectorsCurves.R
│   ├── coreMets.R
│   ├── familyCore.R
│   ├── globalVolcano.R
│   ├── its2.R
│   ├── map.R
│   ├── metsOverview.R
│   ├── pcoa.R
│   └── ml/
├── cotw.Rproj
└── README.md
```

## Analysis scripts

| Script               | Description                                                                                                        |
| -------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `collectorsCurves.R` | Generates metabolite accumulation and collector’s curves                                                          |
| `coreMets.R`         | Identifies and analyzes the convergent core of 66 metabolites                                                     |
| `familyCore.R`       | Identifies and analyzes metabolites detected across all sampled Scleractinian families                            |
| `globalVolcano.R`    | Tests and visualizes global enrichment patterns across compound classes                                           |
| `its2.R`             | Analyzes relationships between ITS2-defined symbiont genera and coral metabolomic composition                |
| `map.R`              | Produces map of sampling locations                                                                               |
| `metsOverview.R`     | Summarizes metabolite annotations and compound-class assignments                                                  |
| `pcoa.R`             | Performs PCoA, PERMANOVA, dendrograms of metabolomic relationships             |
| `ml/`                | Machine learning model files and scripts (see below) |

## Machine learning analyses  

The `analysis/ml/` directory contains:
* Serialized `.joblib` files for trained random forest and gradient boosting models  
* `hpt.py` code for hyperparameter tuning and training random forest and gradient boosting models  
* `inference.ipynb` code for evaluating the trained models on test set and outputting impurity-based feature importance metrics  
* `importance.R` code for visualizing important metabolite features 

## Contact

Contact hs325@duke.edu or smokinroachjr@gmail.com with questions. Citation information will be added upon study publication.  

