#!/usr/bin/perl

## truth
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne_inf.txt -t genarch_known.txt -o so_truth.txt -n 1 -s 0\n";

## drift only 
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne.txt -t genarch_known.txt -o so_drift.txt -n 100 -s 0\n";

## uncertainty in selection (weather known)
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wknown_$i.txt -f ne_inf.txt -t genarch_known.txt -c 0.1 -d 0.1 -o tso_selection$i.txt -n 1 -s 0\n";
}
system "cat tso_selection* > so_selection.txt\n";

## uncertainty in weather + selection 
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne_inf.txt -t genarch_known.txt -c 0.1 -d 0.1 -o tso_weather$i.txt -n 1 -s 0\n";
}
system "cat tso_weather* > so_weather.txt\n";

## uncertainty in weather + selection + drift
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne.txt -t genarch_known.txt -c 0.1 -d 0.1 -o tso_dr_weather$i.txt -n 1 -s 0\n";
}
system "cat tso_dr_weather* > so_dr_weather.txt\n";

## gen uncertain

## drift only 
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne.txt -t genarch_unc.txt -o so_gunc_drift.txt -n 100 -s 0\n";

## uncertainty in selection (weather known)
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wknown_$i.txt -f ne_inf.txt -t genarch_unc.txt -c 0.1 -d 0.1 -o tso_gunc_selection$i.txt -n 1 -s 0\n";
}
system "cat tso_gunc_selection* > so_gunc_selection.txt\n";

## uncertainty in weather + selection 
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne_inf.txt -t genarch_unc.txt -c 0.1 -d 0.1 -o tso_gunc_weather$i.txt -n 1 -s 0\n";
}
system "cat tso_gunc_weather* > so_gunc_weather.txt\n";

## uncertainty in weather + selection + drift
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne.txt -t genarch_unc.txt -c 0.1 -d 0.1 -o tso_gunc_dr_weather$i.txt -n 1 -s 0\n";
}
system "cat tso_gunc_dr_weather* > so_gunc_dr_weather.txt\n";
