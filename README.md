# VibrioCARE

"Signature genomic traits of the core and accessory genome of Vibrio Cholerae O1 drive lineage transmission and disease severity" by Alexandre Maciel-Guerra, Kubra Babaarslan, Michelle Baker, Aura Rahman, Maqsud Muhammad Hossain, Abdus Sadique, Jahidul Alam, Salim Uzzaman, Mohammad Ferdous Rahman Sarker, Nasrin Sultana, Ashraful Islam Khan, Yasmin Ara Begum, Mokibul Hassan Afrad, Nicola Senin, Zakir Hossain Habib, Tahmina Shirin, Firdausi Qadri, and Tania Dottorini

Any questions should be made to the corresponding author Dr Tania Dottorini (Tania.Dottorini@nottingham.ac.uk)

Four scripts are available:

1. Fisher exact test for the lineage separation and Mann Whitney U tests for Figure S4: Lineage_Separation.py
2. SNP network (Figure 2): snp_network.py
3. The classification framework for the binary clinical symptoms (Vomit and Abdominal Pain): ML_pipeline_binary.py
4. The classification framework for the multiclass clinical symptoms (Duration of diarrhoea, Number of stools in 24 hours, Dehydration) - a One-vs-One approach was used to change the multiclass problem to a binary problem: ML_pipeline_multi.py

# System Requirements

## Software requirements

### OS Requirements

This package is supported for Windows. The package has been tested on the following system: 

* Windows: Windows 11 Pro version 23H2 OS build 22631.3296 Windows Feature Experience Pack 1000.22687.1000.0 


### Python Dependencies

```
python v3.9.15
numpy v1.21.5
pandas v1.4.4
scikit-learn v1.2.1
scipy v1.9.3
networkx v2.8.4
matplotlib v3.6.2
```

# Installation Guide:

## Install from Github

```
git clone https://github.com/tan0101/VibrioCARE
cd VibrioCARE
python3 setup.py install
```

# Algorithm's Flow

# License

This project is covered under the **AGPL-3.0 license**.
