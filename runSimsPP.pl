#!/usr/bin/perl

## truth
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne_inf.txt -t genarch_known.txt -o sop_truth.txt -n 1 -s 1\n";

## drift only 
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne.txt -t genarch_known.txt -o sop_drift.txt -n 100 -s 1\n";

## uncertainty in selection (weather known)
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wknown_$i.txt -f ne_inf.txt -t genarch_known.txt -c 0.1 -d 0.1 -o tsop_selection$i.txt -n 1 -s 1\n";
}
system "cat tsop_selection* > sop_selection.txt\n";

## uncertainty in weather + selection 
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne_inf.txt -t genarch_known.txt -c 0.1 -d 0.1 -o tsop_weather$i.txt -n 1 -s 1\n";
}
system "cat tsop_weather* > sop_weather.txt\n";

## uncertainty in weather + selection + drift
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne.txt -t genarch_known.txt -c 0.1 -d 0.1 -o tsop_dr_weather$i.txt -n 1 -s 1\n";
}
system "cat tsop_dr_weather* > sop_dr_weather.txt\n";

## gen uncertain

## drift only 
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne.txt -t genarch_unc.txt -o sop_gunc_drift.txt -n 100 -s 1\n";

## uncertainty in selection (weather known)
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wknown_$i.txt -f ne_inf.txt -t genarch_unc.txt -c 0.1 -d 0.1 -o tsop_gunc_selection$i.txt -n 1 -s 1\n";
}
system "cat tsop_gunc_selection* > sop_gunc_selection.txt\n";

## uncertainty in weather + selection 
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne_inf.txt -t genarch_unc.txt -c 0.1 -d 0.1 -o tsop_gunc_weather$i.txt -n 1 -s 1\n";
}
system "cat tsop_gunc_weather* > sop_gunc_weather.txt\n";

## uncertainty in weather + selection + drift
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne.txt -t genarch_unc.txt -c 0.1 -d 0.1 -o tsop_gunc_dr_weather$i.txt -n 1 -s 1\n";
}
system "cat tsop_gunc_dr_weather* > sop_gunc_dr_weather.txt\n";
