# VibrioCARE

"Signature genomic traits of the core and accessory genome of Vibrio Cholerae O1 drive lineage transmission and disease severity" by Alexandre Maciel-Guerra, Kubra Babaarslan, Michelle Baker, Aura Rahman, Maqsud Muhammad Hossain, Abdus Sadique, Jahidul Alam, Salim Uzzaman, Mohammad Ferdous Rahman Sarker, Nasrin Sultana, Ashraful Islam Khan, Yasmin Ara Begum, Mokibul Hassan Afrad, Nicola Senin, Zakir Hossain Habib, Tahmina Shirin, Firdausi Qadri, and Tania Dottorini

Any questions should be made to the corresponding author Dr Tania Dottorini (Tania.Dottorini@nottingham.ac.uk)

Four scripts are available:

1. Fisher exact test for the lineage separation and Mann Whitney U tests for Figure S4: Lineage_Separation.py
2. SNP network (Figure 2): snp_network.py
3. The classification framework for the binary clinical symptoms (Vomit and Abdominal Pain): ML_pipeline_binary.py
4. The classification framework for the multiclass clinical symptoms (Duration of diarrhoea, Number of stools in 24 hours, Dehydration) - a One-vs-One approach was used to change the multiclass problem to a binary problem: ML_pipeline_multi.py
