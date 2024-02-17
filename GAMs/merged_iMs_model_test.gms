$title No Activator, metaboltie binding inhibitor
Option Seed=42;

* Load in cEffector Matrix
* Set Dimensions
Set
   iM 'iModulons'
   gene 'genes'
   sample 'samples'
   stable_constant 'constant across all'
       / kd_act_metab,
         kd_inh_metab /
   gene_constant 'grid basal constants'
       / KdRNAP,
         KdRNAPCrp,
         mRNA_total,
         cell_volume,
         k_d_TF,
         kDeg,
         promoterConcVal,
         TF,
         u,
         kEscape,
         KeqOpening /
    sample_constant 'grid basal constants'
       / RNAP /
    limits 'input parameters'
       / act_TF_conc_lo,
         act_TF_conc_up,
         act_Kd_lo,
         act_Kd_up,
         inh_TF_conc_lo,
         inh_TF_conc_up,
         inh_Kd_lo,
         inh_Kd_up,
         weight_act_obj1,
         weight_inh_obj1,
         weight_act_obj2,
         weight_inh_obj2,
         weight_mRNA_match,
         weight_act_corr,
         weight_inh_corr,
         act_metab_Total_lo,
         act_metab_Total_up,
         inh_metab_Total_lo,
         inh_metab_Total_up /

* Load in saved values
$onUNDF
$call csv2gdx ./input_files/composite_cAct_vals.csv id=cAct index=1,2 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_cAct.gdx
$gdxIn ./input_GDX/input_cAct.gdx
Parameter cAct(gene, iM, sample) 'cAct values';
$load cAct
$gdxIn

$call csv2gdx ./input_files/composite_cInh_vals.csv id=cInh index=1,2 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_cInh.gdx
$gdxIn ./input_GDX/input_cInh.gdx
Parameter cInh(gene, iM, sample) 'cInh values';
$load cInh
$gdxIn


Display cInh;
Display gene;
Display iM;
Display sample;