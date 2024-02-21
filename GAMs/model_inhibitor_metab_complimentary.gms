$title No Activator, metaboltie binding inhibitor
Option Seed=42;

* Load in cEffector Matrix
* Set Dimensions
Set
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
         inh_metab_Total_lo,
         inh_metab_Total_up /

* Load in saved values
$call csv2gdx ./input_files/composite_cAct_vals.csv id=cAct index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_cAct.gdx
$gdxIn ./input_GDX/input_cAct.gdx
$load sample = dim1
$load gene = dim2
Parameter cAct(sample, gene) 'cAct values';
$load cAct
$gdxIn

$call csv2gdx ./input_files/composite_cInh_vals.csv id=cInh index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_cInh.gdx
$gdxIn ./input_GDX/input_cInh.gdx
Parameter cInh(sample, gene) 'cInh values';
$load cInh
$gdxIn

$call csv2gdx ./input_files/exported_act_TF_conc.csv id=meas_act_TF index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_act_TF_conc.gdx
$gdxIn ./input_GDX/input_act_TF_conc.gdx
Parameter meas_act_TF(sample) 'actual act TF values';
$load meas_act_TF
$gdxIn

$call csv2gdx ./input_files/exported_inh_TF_conc.csv id=meas_inh_TF index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_inh_TF_conc.gdx
$gdxIn ./input_GDX/input_inh_TF_conc.gdx
Parameter meas_inh_TF(sample) 'actual inh TF values';
$load meas_inh_TF
$gdxIn

$call csv2gdx ./input_files/grid_constants.csv id=basal_constants index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/basal.gdx
$gdxIn ./input_GDX/basal.gdx
Parameter basal_constants(gene_constant, gene) 'basal constants';
$load basal_constants
$gdxIn

$call csv2gdx ./input_files/sample_constants.csv id=sample_constants index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/sample_basal.gdx
$gdxIn ./input_GDX/sample_basal.gdx
Parameter sample_constants(sample_constant, sample) 'sample constants';
$load sample_constants
$gdxIn

$call csv2gdx ./input_files/input_constants.csv id=input_constants index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_constants.gdx
$gdxIn ./input_GDX/input_constants.gdx
Parameter input_constants(stable_constant) 'input constants';
$load input_constants
$gdxIn

$call csv2gdx ./input_files/actual_mRNA_ratio.csv id=actual_mRNA index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/actual_mRNA.gdx
$gdxIn ./input_GDX/actual_mRNA.gdx
Parameter actual_mRNA(sample, gene) 'actual_mRNA';
$load actual_mRNA
$gdxIn

$call csv2gdx ./input_files/max_cActs.csv id=max_cAct index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/max_cAct.gdx
$gdxIn ./input_GDX/max_cAct.gdx
Parameter max_cAct(gene) 'max_cAct';
$load max_cAct
$gdxIn

$call csv2gdx ./input_files/max_cInhs.csv id=max_cInh index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/max_cInh.gdx
$gdxIn ./input_GDX/max_cInh.gdx
Parameter max_cInh(gene) 'max_cInh';
$load max_cInh
$gdxIn

$call csv2gdx ./input_files/parameters.csv id=input_parameters index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_parameters.gdx
$gdxIn ./input_GDX/input_parameters.gdx
Parameter input_parameters(limits) 'input_paras';
$load input_parameters
$gdxIn

***
* Set up other matrices
Variables
    act_TF_conc(sample) 'TF concentration'
    act_Kd(gene) 'Kd Values'
    inh_Kd(gene) 'Kd Values'
    act_metab_Total(sample)
    inh_metab_Total(sample);

Display input_parameters;
Display sample_constants;

* constrain matrices
act_TF_conc.lo(sample) = log10(input_parameters('act_TF_conc_lo'));
act_TF_conc.up(sample) = log10(input_parameters('act_TF_conc_up'));
act_Kd.lo(gene) = log10(input_parameters('act_Kd_lo'));
act_Kd.up(gene) = log10(input_parameters('act_Kd_up'));
inh_Kd.lo(gene) = log10(input_parameters('inh_Kd_lo'));
inh_Kd.up(gene) = log10(input_parameters('inh_Kd_up'));
inh_metab_Total.lo(sample) = log10(input_parameters('inh_metab_Total_lo'));
inh_metab_Total.up(sample) = log10(input_parameters('inh_metab_Total_up'));


$onText
* standard values for above
act_TF_conc.lo(sample) = log10(1E-8);
act_TF_conc.up(sample) = log10(1E-5);
act_Kd.lo(gene) = log10(1E-10);
act_Kd.up(gene) = log10(1E-6);
inh_Kd.lo(gene) = log10(1E-10);
inh_Kd.up(gene) = log10(1E-6);
act_metab_Total.lo(gene) = log10(1E-10);
act_metab_Total.up(gene) = log10(1E-6);
inh_metab_Total.lo(gene) = log10(1E-5);
inh_metab_Total.up(gene) = log10(1E-3);
$offText

* initialize matrices
act_TF_conc.l(sample) = uniformInt(-10, -3);
act_Kd.l(gene) = uniformInt(-10, -3);
inh_Kd.l(gene) = uniformInt(-10, -3);
inh_metab_Total.l(sample) = uniformInt(-10, -3);

* Now set up the optimization
Variables
    act_diff1   difference between TF_conc div Kd and Ceff
    inh_diff1   similar
    match_diff
    total_diff;

* calculated sub variables
Variables
    cAct_calc(sample, gene)
    cInh_calc(sample, gene);

Equations
    eq_cAct_calc(sample, gene)
    eq_cInh_calc(sample, gene);

Equations
    act_obj1   difference between TF_conc div Kd and cEff
    inh_obj1   difference between TF_conc div Kd and cEff
    match_obj  difference between predicted and actual mRNA
    total_obj;


* Set weights for the two objectives (you can adjust these weights as needed)
* mRNA weighting
Scalar weight_balance3 /195260.5676/;
* TF conenctration correlation matching
Scalar weight_balance4 /1953.70798879785/;
Scalar weight_act_obj1 /1/;
Scalar weight_inh_obj1 /1/;
Scalar weight_mRNA_match /0/;
Scalar weight_act_corr /0/;
Scalar weight_inh_corr /0/;


* set them from parameters file
weight_act_obj1 = input_parameters('weight_act_obj1');
weight_inh_obj1 = input_parameters('weight_inh_obj1');
weight_mRNA_match = input_parameters('weight_mRNA_match');
weight_act_corr = input_parameters('weight_act_corr');
weight_inh_corr = input_parameters('weight_inh_corr');


* create equations
total_obj .. total_diff =e= weight_balance3 * weight_mRNA_match * match_diff + weight_act_obj1 * act_diff1 + weight_inh_obj1 * inh_diff1;


* equations for cInhibitor and cActivator (cActivator is basically null and unused right now)
eq_cAct_calc(sample, gene) .. cAct_calc(sample, gene) =e= 10**act_TF_conc(sample) / 10**act_Kd(gene);

eq_cInh_calc(sample, gene) .. cInh_calc(sample, gene) =e= (
        3 * 10**inh_metab_Total(sample) * 10**inh_Kd(gene) + 
        input_constants('kd_inh_metab') * 10**inh_Kd(gene) + 
        3 * 10**inh_Kd(gene) * meas_inh_TF(sample) + 
        (
            -36 * 10**inh_metab_Total(sample) * (10**inh_Kd(gene))**2 * meas_inh_TF(sample) + 
             (
                 3 * 10**inh_metab_Total(sample) * 10**inh_Kd(gene) + 
                  input_constants('kd_inh_metab') * 10**inh_Kd(gene) + 
                  3 * 10**inh_Kd(gene) * meas_inh_TF(sample)
             )**2
        )**(.5)
    )
    / (18 * (10**inh_Kd(gene))**2);
    
* objective equations
act_obj1 .. act_diff1 =e= sum((gene, sample), (abs((cAct_calc(sample, gene) - cAct(sample, gene))) / max_cAct(gene) )**2);

inh_obj1 .. inh_diff1 =e= sum((gene, sample), (abs((cInh_calc(sample, gene) - cInh(sample, gene))) / max_cInh(gene) )**2);


match_obj .. match_diff =e= sum((gene, sample), (abs(actual_mRNA(sample, gene) - ((cAct_calc(sample, gene)*basal_constants('KdRNAP', gene) + basal_constants('KdRNAPCrp', gene))*(basal_constants('KdRNAP', gene) + sample_constants('RNAP', sample) +  basal_constants('KeqOpening', gene)*sample_constants('RNAP', sample))) / (((1 + cAct_calc(sample, gene) + cInh_calc(sample, gene))*basal_constants('KdRNAP', gene)*basal_constants('KdRNAPCrp', gene) + cAct_calc(sample, gene)*basal_constants('KdRNAP', gene)*(1 + basal_constants('KeqOpening', gene))*sample_constants('RNAP', sample) + basal_constants('KdRNAPCrp', gene)*(1 + basal_constants('KeqOpening', gene))*sample_constants('RNAP', sample)))))**2 );



* modify model parameters
Option Iterlim=10000;
*$set dnlp acc=1e-15  // Set the accuracy or tolerance level
*$set dnlp step=1e-15  // Set the step size

* run the model
Model ElementWiseOptimization /all/;
ElementWiseOptimization.optfile = 1;
Solve ElementWiseOptimization using dnlp minimizing total_diff;


* Export results
execute_unload "./output_GDX/cAct_TF_conc_results.gdx" act_TF_conc.L act_TF_conc.M
execute 'gdxdump ./output_GDX/cAct_TF_conc_results.gdx noData > ./output_files/cAct_TF_conc_results.csv symb=act_TF_conc format=csv';

execute_unload "./output_GDX/cAct_Kd_results.gdx" act_Kd.L act_Kd.M
execute 'gdxdump ./output_GDX/cAct_Kd_results.gdx noData > ./output_files/cAct_Kd_results.csv symb=act_Kd format=csv';

execute_unload "./output_GDX/cInh_Kd_results.gdx" inh_Kd.L inh_Kd.M
execute 'gdxdump ./output_GDX/cInh_Kd_results.gdx noData > ./output_files/cInh_Kd_results.csv symb=inh_Kd format=csv';

execute_unload "./output_GDX/inh_metab_Total.gdx" inh_metab_Total.L inh_metab_Total.M
execute 'gdxdump ./output_GDX/inh_metab_Total.gdx noData > ./output_files/inh_metab_Total.csv symb=inh_metab_Total format=csv';