#!/usr/bin/perl

## truth
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne_inf.txt -t genarch_known.txt -o soa_truth.txt -n 1 -s 0\n";

## drift only 
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne.txt -t genarch_known.txt -o soa_drift.txt -n 100 -s 0\n";

## uncertainty in selection (weather known)
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wknown_$i.txt -f ne_inf.txt -t genarch_known.txt -c 0.058 -d 0.058 -o tsoa_selection$i.txt -n 1 -s 0\n";
}
system "cat tsoa_selection* > soa_selection.txt\n";

## uncertainty in weather + selection 
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne_inf.txt -t genarch_known.txt -c 0.058 -d 0.058 -o tsoa_weather$i.txt -n 1 -s 0\n";
}
system "cat tsoa_weather* > soa_weather.txt\n";

## uncertainty in weather + selection + drift
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne.txt -t genarch_known.txt -c 0.058 -d 0.058 -o tsoa_dr_weather$i.txt -n 1 -s 0\n";
}
system "cat tsoa_dr_weather* > soa_dr_weather.txt\n";

## gen uncertain

## drift only 
system "./psims -g p0.txt -e envres_wknown_1.txt -f ne.txt -t genarch_unc.txt -o soa_gunc_drift.txt -n 100 -s 0\n";

## uncertainty in selection (weather known)
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wknown_$i.txt -f ne_inf.txt -t genarch_unc.txt -c 0.058 -d 0.058 -o tsoa_gunc_selection$i.txt -n 1 -s 0\n";
}
system "cat tsoa_gunc_selection* > soa_gunc_selection.txt\n";

## uncertainty in weather + selection 
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne_inf.txt -t genarch_unc.txt -c 0.058 -d 0.058 -o tsoa_gunc_weather$i.txt -n 1 -s 0\n";
}
system "cat tsoa_gunc_weather* > soa_gunc_weather.txt\n";

## uncertainty in weather + selection + drift
foreach $i (1..100){
	system "sleep 2\n";
	system "./psims -g p0.txt -e envres_wunc_$i.txt -f ne.txt -t genarch_unc.txt -c 0.058 -d 0.058 -o tsoa_gunc_dr_weather$i.txt -n 1 -s 0\n";
}
system "cat tsoa_gunc_dr_weather* > soa_gunc_dr_weather.txt\n";
