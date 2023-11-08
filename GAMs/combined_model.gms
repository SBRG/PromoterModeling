$title cActivator toy example

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
         KeqOpening /;

* Load in saved values
$call csv2gdx ../data/save_for_GAMs/composite_cAct_vals.csv id=cAct index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_cAct.gdx
$gdxIn ./input_GDX/input_cAct.gdx
$load sample = dim1
$load gene = dim2
Parameter cAct(sample, gene) 'cAct values';
$load cAct
$gdxIn

$call csv2gdx ../data/save_for_GAMs/composite_cInh_vals.csv id=cInh index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_cInh.gdx
$gdxIn ./input_GDX/input_cInh.gdx
Parameter cInh(sample, gene) 'cInh values';
$load cInh
$gdxIn

$call csv2gdx ../data/save_for_GAMs/exported_act_TF_conc.csv id=meas_act_TF index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_act_TF_conc.gdx
$gdxIn ./input_GDX/input_act_TF_conc.gdx
Parameter meas_act_TF(sample) 'actual act TF values';
$load meas_act_TF
$gdxIn

$call csv2gdx ../data/save_for_GAMs/exported_inh_TF_conc.csv id=meas_inh_TF index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/input_inh_TF_conc.gdx
$gdxIn ./input_GDX/input_inh_TF_conc.gdx
Parameter meas_inh_TF(sample) 'actual inh TF values';
$load meas_inh_TF
$gdxIn

$call csv2gdx ../data/save_for_GAMs/grid_constants.csv id=basal_constants index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/basal.gdx
$gdxIn ./input_GDX/basal.gdx
Parameter basal_constants(constant, gene) 'basal constants';
$load basal_constants
$gdxIn

$call csv2gdx ../data/save_for_GAMs/actual_mRNA_ratio.csv id=actual_mRNA index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/actual_mRNA.gdx
$gdxIn ./input_GDX/actual_mRNA.gdx
Parameter actual_mRNA(sample, gene) 'actual_mRNA';
$load actual_mRNA
$gdxIn

$call csv2gdx ../data/save_for_GAMs/max_cActs.csv id=max_cAct index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/max_cAct.gdx
$gdxIn ./input_GDX/max_cAct.gdx
Parameter max_cAct(gene) 'max_cAct';
$load max_cAct
$gdxIn

$call csv2gdx ../data/save_for_GAMs/max_cInhs.csv id=max_cInh index=1 values=2..lastCol useHeader=y trace=0 output=./input_GDX/max_cInh.gdx
$gdxIn ./input_GDX/max_cInh.gdx
Parameter max_cInh(gene) 'max_cInh';
$load max_cInh
$gdxIn

***
* Set up other matrices
Variables
    act_TF_conc(sample) 'TF concentration'
    act_Kd(gene) 'Kd Values'
    inh_TF_conc(sample) 'TF concentration'
    inh_Kd(gene) 'Kd Values';

* constrain matrices
act_TF_conc.lo(sample) = log10(1E-15);
act_TF_conc.up(sample) = 5+log10(meas_act_TF(sample));
act_Kd.lo(gene) = log10(1E-9);
act_Kd.up(gene) = log10(1E-4);
inh_TF_conc.lo(sample) = log10(1E-15);
inh_TF_conc.up(sample) = 5+log10(meas_inh_TF(sample));
inh_Kd.lo(gene) = log10(1E-9);
inh_Kd.up(gene) = log10(1E-4);

* initialize matrices
act_TF_conc.l(sample) = uniformInt(1, 10);
act_Kd.l(gene) = uniformInt(1, 10);
inh_TF_conc.l(sample) = uniformInt(1, 10);
inh_Kd.l(gene) = uniformInt(1, 10);

* Now set up the optimization
Variables
    act_diff1   difference between TF_conc div Kd and Ceff
    inh_diff1   similar
    act_diff2   difference between measured TF conc and predicted TF conc
    inh_diff2   similar
    match_diff
    total_diff;

Equations
    act_obj1   difference between TF_conc div Kd and cEff
    inh_obj1   difference between TF_conc div Kd and cEff
    act_obj2   difference between measured TF conc and predicted TF conc
    inh_obj2   difference between measured TF conc and predicted TF conc
    match_obj  difference between predicted and actual mRNA
    total_obj;

* Set weights for the two objectives (you can adjust these weights as needed)
Scalar weight_balance2 /230632.01555323132/;
Scalar weight_balance3 /0.1874203277944165/;
Scalar weight_act_obj1 /1/;
Scalar weight_inh_obj1 /1/;
Scalar weight_act_obj2 /0/;
Scalar weight_inh_obj2 /0/;
Scalar weight_mRNA_match /.5/;
* ^ obj1 is just way bigger by how its calculated, this rebalances them to be about even and obj1_2 can scale them to be different

total_obj .. total_diff =e= weight_balance3 * weight_mRNA_match * match_diff + weight_act_obj1 * act_diff1 + weight_inh_obj1 * inh_diff1 + weight_balance2 * (weight_act_obj2 * act_diff2 + weight_inh_obj2 * inh_diff2);

act_obj1 .. act_diff1 =e= sum((gene, sample), abs((10**act_TF_conc(sample) / 10**act_Kd(gene) - cAct(sample, gene))/max_cAct(gene)));

inh_obj1 .. inh_diff1 =e= sum((gene, sample), abs((10**inh_TF_conc(sample) / 10**inh_Kd(gene) - cInh(sample, gene))/max_cInh(gene)));

act_obj2 .. act_diff2 =e= sum(sample, (10**(meas_act_TF(sample)) - 0.1*10**act_TF_conc(sample))**2);

inh_obj2 .. inh_diff2 =e= sum(sample, (10**(meas_inh_TF(sample)) - 0.1*10**inh_TF_conc(sample))**2);

match_obj .. match_diff =e= sum((gene, sample), abs(actual_mRNA(sample, gene) - (((10**act_TF_conc(sample) / 10**act_Kd(gene))*basal_constants('KdRNAP', gene) + basal_constants('KdRNAPCrp', gene))*(basal_constants('KdRNAP', gene) + basal_constants('RNAP', gene) +  basal_constants('KeqOpening', gene)*basal_constants('RNAP', gene))) / ((1 + (10**act_TF_conc(sample) / 10**act_Kd(gene)) + (10**inh_TF_conc(sample) / 10**inh_Kd(gene)))*basal_constants('KdRNAP', gene)*basal_constants('KdRNAPCrp', gene) + (10**act_TF_conc(sample) / 10**act_Kd(gene))*basal_constants('KdRNAP', gene)*(1 + basal_constants('KeqOpening', gene))*basal_constants('RNAP', gene) + basal_constants('KdRNAPCrp', gene)*(1 + basal_constants('KeqOpening', gene))*basal_constants('RNAP', gene))))


*need to add existing mRNA value to the above equation


* ^ I found one paper that said about 10% of crp is active based on cAMP presence, this should be changed later though

* modify model parameters
*$set dnlp maxiter=10000  // Increase the maximum number of iterations
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
execute 'gdxdump ./output_GDX/cAct_TF_conc_results.gdx noData > ../data/GAMS_output/cAct_TF_conc_results.csv symb=act_TF_conc format=csv';
execute 'gdxdump ./output_GDX/cAct_Kd_results.gdx noData > ../data/GAMS_output/cAct_Kd_results.csv symb=act_Kd format=csv';
execute_unload "./output_GDX/cInh_TF_conc_results.gdx" inh_TF_conc.L inh_TF_conc.M
execute_unload "./output_GDX/cInh_Kd_results.gdx" inh_Kd.L inh_Kd.M
execute 'gdxdump ./output_GDX/cInh_TF_conc_results.gdx noData > ../data/GAMS_output/cInh_TF_conc_results.csv symb=inh_TF_conc format=csv';
execute 'gdxdump ./output_GDX/cInh_Kd_results.gdx noData > ../data/GAMS_output/cInh_Kd_results.csv symb=inh_Kd format=csv';
