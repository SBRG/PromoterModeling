
$gdxIn all.gdx
$onEmpty
$offEolCom
$eolCom !!

Set gene(*) genes ;
$loadDC gene

Set sample(*) samples ;
$loadDC sample

Parameter cEff(sample,gene) cEff values ;
$loadDC cEff

Parameter meas_TF(sample) actual TF values ;
$loadDC meas_TF

free     Variable TF_conc(sample) TF concentration ;
$loadDC TF_conc

free     Variable Kd(gene) Kd Values ;
$loadDC Kd

free     Variable diff1 difference between TF_conc div Kd and Ceff ;
$loadDC diff1

free     Variable diff2 difference between measured TF conc and predicted TF conc ;
$loadDC diff2

free     Variable total_diff ;
$loadDC total_diff

Equation obj1 difference between TF_conc div Kd and cEff ;
$loadDC obj1

Equation obj2 difference between measured TF conc and predicted TF conc ;
$loadDC obj2

Equation total_obj ;
$loadDC total_obj

Scalar weight_obj1 ;
$loadDC weight_obj1

Scalar weight_obj2 ;
$loadDC weight_obj2

Set days(*) ;
$loadDC days

Set weekend(days) ;
$loadDC weekend

Scalar DAYSINWEEK ;
$loadDC DAYSINWEEK

Set item(*) ;
$loadDC item

Parameter ITEMCOST(item) ;
$loadDC ITEMCOST

$offEmpty
