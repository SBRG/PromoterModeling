$title scaling the combined_model to be able to run more than one promoter at a time

* Load in cEffector Matrix
* Set Dimensions
Set
   gene 'genes'
   sample 'samples'
   constant 'grid basal constants'
       / KdRNAP,
         KdRNAPCrp,
         RNAP,
         mRNA_total,
         cell_volume,
         k_d_TF,
         kDeg,
         promoterConcVal,
         TF,
         u,
         kEscape,
         KeqOpening /
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
         weight_inh_corr /

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
Parameter basal_constants(constant, gene) 'basal constants';
$load basal_constants
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
    inh_TF_conc(sample) 'TF concentration'
    inh_Kd(gene) 'Kd Values';

* constrain matrices
act_TF_conc.lo(sample) = log10(input_parameters('act_TF_conc_lo'));
act_TF_conc.up(sample) = log10(input_parameters('act_TF_conc_up'));
*log10(meas_act_TF(sample));
act_Kd.lo(gene) = log10(input_parameters('act_Kd_lo'));
act_Kd.up(gene) = log10(input_parameters('act_Kd_up'));
inh_TF_conc.lo(sample) = log10(input_parameters('inh_TF_conc_lo'));
inh_TF_conc.up(sample) = log10(input_parameters('inh_TF_conc_up'));
*log10(meas_inh_TF(sample));
inh_Kd.lo(gene) = log10(input_parameters('inh_Kd_lo'));
inh_Kd.up(gene) = log10(input_parameters('inh_Kd_up'));
$onText
* standard values for above
act_TF_conc.lo(sample) = log10(1E-8);
act_TF_conc.up(sample) = log10(1E-5);
act_Kd.lo(gene) = log10(1E-10);
act_Kd.up(gene) = log10(1E-6);
inh_TF_conc.lo(sample) = log10(1E-8);
inh_TF_conc.up(sample) = log10(1E-5);
inh_Kd.lo(gene) = log10(1E-10);
inh_Kd.up(gene) = log10(1E-6);
$offText

* initialize matrices
act_TF_conc.l(sample) = uniformInt(-10, -3);
act_Kd.l(gene) = uniformInt(-10, -3);
inh_TF_conc.l(sample) = uniformInt(-10, -3);
inh_Kd.l(gene) = uniformInt(-10, -3);

* Now set up the optimization
Variables
    act_diff1   difference between TF_conc div Kd and Ceff
    inh_diff1   similar
    act_diff2   difference between measured TF conc and predicted TF conc
    inh_diff2   similar
    match_diff
    act_tf_conc_correlation_diff
    inh_tf_conc_correlation_diff
    act_avg_actual_TF_conc
    act_avg_pred_TF_conc
    inh_avg_actual_TF_conc
    inh_avg_pred_TF_conc
    total_diff;

Equations
    act_avg_actual_TF_conc_eq
    act_avg_pred_TF_conc_eq
    inh_avg_actual_TF_conc_eq
    inh_avg_pred_TF_conc_eq;
    
Equations
    act_obj1   difference between TF_conc div Kd and cEff
    inh_obj1   difference between TF_conc div Kd and cEff
    act_obj2   difference between measured TF conc and predicted TF conc
    inh_obj2   difference between measured TF conc and predicted TF conc
    match_obj  difference between predicted and actual mRNA
    act_tf_conc_correlation_obj
    inh_tf_conc_correlation_obj
    total_obj;

* Set weights for the two objectives (you can adjust these weights as needed)
* TF Concentration match weighting
Scalar weight_balance2 /9422806.40775587/;
*mRNA weighting
Scalar weight_balance3 /60.56763055777653/;
* TF conenctration correlation matching
Scalar weight_balance4 /1953.70798879785/;
Scalar weight_act_obj1 /1/;
Scalar weight_inh_obj1 /1/;
Scalar weight_act_obj2 /0/;
Scalar weight_inh_obj2 /0/;
Scalar weight_mRNA_match /0/;
Scalar weight_act_corr /0/;
Scalar weight_inh_corr /0/;

* set them from parameters file
weight_act_obj1 = input_parameters('weight_act_obj1');
weight_inh_obj1 = input_parameters('weight_inh_obj1');
weight_act_obj2 = input_parameters('weight_act_obj2');
weight_inh_obj2 = input_parameters('weight_inh_obj2');
weight_mRNA_match = input_parameters('weight_mRNA_match');
weight_act_corr = input_parameters('weight_act_corr');
weight_inh_corr = input_parameters('weight_inh_corr');


* create equations
total_obj .. total_diff =e= weight_balance4 * (weight_act_corr * act_tf_conc_correlation_diff + weight_inh_corr * inh_tf_conc_correlation_diff) + weight_balance3 * weight_mRNA_match * match_diff + weight_act_obj1 * act_diff1 + weight_inh_obj1 * inh_diff1 + weight_balance2 * (weight_act_obj2 * act_diff2 + weight_inh_obj2 * inh_diff2);

act_obj1 .. act_diff1 =e= sum((gene, sample), abs((10**act_TF_conc(sample) / 10**act_Kd(gene) - cAct(sample, gene))));

inh_obj1 .. inh_diff1 =e= sum((gene, sample), abs((10**inh_TF_conc(sample) / 10**inh_Kd(gene) - cInh(sample, gene))));

act_obj2 .. act_diff2 =e= sum(sample, abs(meas_act_TF(sample) - 0.1*(10**act_TF_conc(sample)))**2);

inh_obj2 .. inh_diff2 =e= sum(sample, abs(meas_inh_TF(sample) - 0.1*(10**inh_TF_conc(sample)))**2);

match_obj .. match_diff =e= sum((gene, sample), abs(actual_mRNA(sample, gene) - (((10**act_TF_conc(sample) / 10**act_Kd(gene))*basal_constants('KdRNAP', gene) + basal_constants('KdRNAPCrp', gene))*(basal_constants('KdRNAP', gene) + basal_constants('RNAP', gene) +  basal_constants('KeqOpening', gene)*basal_constants('RNAP', gene))) / ((1 + (10**act_TF_conc(sample) / 10**act_Kd(gene)) + (10**inh_TF_conc(sample) / 10**inh_Kd(gene)))*basal_constants('KdRNAP', gene)*basal_constants('KdRNAPCrp', gene) + (10**act_TF_conc(sample) / 10**act_Kd(gene))*basal_constants('KdRNAP', gene)*(1 + basal_constants('KeqOpening', gene))*basal_constants('RNAP', gene) + basal_constants('KdRNAPCrp', gene)*(1 + basal_constants('KeqOpening', gene))*basal_constants('RNAP', gene))));

* This is Pearson correlation
act_avg_actual_TF_conc_eq .. act_avg_actual_TF_conc =e= sum(sample, meas_act_TF(sample)) / sum(sample, 1);
act_avg_pred_TF_conc_eq .. act_avg_pred_TF_conc =e= sum(sample, 10**act_TF_conc(sample)) / sum(sample, 1);
inh_avg_actual_TF_conc_eq .. inh_avg_actual_TF_conc =e= sum(sample, meas_inh_TF(sample)) / sum(sample, 1);
inh_avg_pred_TF_conc_eq .. inh_avg_pred_TF_conc =e= sum(sample, 10**inh_TF_conc(sample)) / sum(sample, 1);

act_tf_conc_correlation_obj .. act_tf_conc_correlation_diff =e= -1*(sum(sample, (10**act_TF_conc(sample) - act_avg_pred_TF_conc)*(meas_act_TF(sample) - act_avg_actual_TF_conc)) / ((sum(sample, abs(10**act_TF_conc(sample) - act_avg_pred_TF_conc)**2))**.5 * (sum(sample, abs(meas_act_TF(sample) - act_avg_actual_TF_conc)**2))**.5));

inh_tf_conc_correlation_obj .. inh_tf_conc_correlation_diff =e= sum(sample, (10**inh_TF_conc(sample) - inh_avg_pred_TF_conc)*(meas_inh_TF(sample) - inh_avg_actual_TF_conc)) / ((sum(sample, abs(10**inh_TF_conc(sample) - inh_avg_pred_TF_conc)**2))**.5 * (sum(sample, abs(meas_inh_TF(sample) - inh_avg_actual_TF_conc)**2))**.5);



* ^ I found one paper that said about 10% of crp is active based on cAMP presence, this should be changed later though

* modify model parameters
Option Iterlim=10000;
*$set dnlp acc=1e-15  // Set the accuracy or tolerance level
*$set dnlp step=1e-15  // Set the step size

* run the model
Model ElementWiseOptimization /all/;
Solve ElementWiseOptimization using dnlp minimizing total_diff;

* display results
*Display TF_conc.l, TF_conc.m;
*Display diff1.l, diff2.l;


* Export results
execute_unload "./output_GDX/cAct_TF_conc_results.gdx" act_TF_conc.L act_TF_conc.M
execute_unload "./output_GDX/cAct_Kd_results.gdx" act_Kd.L act_Kd.M
execute 'gdxdump ./output_GDX/cAct_TF_conc_results.gdx noData > ./output_files/cAct_TF_conc_results.csv symb=act_TF_conc format=csv';
execute 'gdxdump ./output_GDX/cAct_Kd_results.gdx noData > ./output_files/cAct_Kd_results.csv symb=act_Kd format=csv';
execute_unload "./output_GDX/cInh_TF_conc_results.gdx" inh_TF_conc.L inh_TF_conc.M
execute_unload "./output_GDX/cInh_Kd_results.gdx" inh_Kd.L inh_Kd.M
execute 'gdxdump ./output_GDX/cInh_TF_conc_results.gdx noData > ./output_files/cInh_TF_conc_results.csv symb=inh_TF_conc format=csv';
execute 'gdxdump ./output_GDX/cInh_Kd_results.gdx noData > ./output_files/cInh_Kd_results.csv symb=inh_Kd format=csv';
