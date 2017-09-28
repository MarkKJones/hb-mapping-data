HB mapping data
================

A variety of mapping data was taken for HB by Eric Sun, Steve Lassiter 
and Mike Fowler from roughly Dec 2015-Mar 2016.

The HB at this time did not have the front clamp and extra side iron which was added later to reduce the
stray fields.

TOSCA calculations were taken from the spreadsheets.

On Sept 27, 2017, Steve Lassiter did TOSCA calculations with model HB-AB which is the configuration at the time
of the Dec 2015-Mar 2016. Steve also did TOSCA calculations with model HB-V10 which is the configuration with addition iron
and front field clamp.

## Hall Probe info

  * A Lakeshore probe Lakeshore-XHGT-9060-55937-Q1 was used for the measurements. This probe is now in Q1.
  The present probe in HB is Lakeshore-XHGT-9060-55587-strayfield-presently-in-HB.pdf which
  was also used for the stray field measurements.
  * The probe correction table is in Lakeshore-XHGT-9060-55937-Q1.pdf .
  * The data from the table has been put in file: Lakeshore-probe.dat
  * This linearity error is different for positive and negative field directions. 
  * The general Lakeshore manual is file: Lakeshore_general_manual.pdf
  * In section 4.5 in the Lakeshore manual it states the the error in a magnitude error.
  1. In the data table, the Error = abs(Measured Field)-abs(True Field). 
  2. So for positive fields True Field = Measured Field - Error. 
  3. For negative fields True Field = Measured Field + Error.

## Data taking description

The field mapping data taking is described in three powerpoint files.

1. The file HB Field Mapping.ppt describes the measurements of the field
along the center line of the HB and stray field measurement along the beam line.
The rig (Rig 1) is shown in the  figure on slide 3 of the ppt file. 

2. The file HB Field Mapping-V3.ppt describes data taken with a different rig
so that measurements could be made along Z at five different X,Y positions.
The rig  (Rig 2) is shown in the  figure on slide 3 of the ppt file. 

3. The file Harmonic Analysis of SHMS HB.ppt describes data taken with another rig
so that measurements could be at center Z at points around a 5.908cm circle.
The rig  (Rig 3) is shown in the  figure on slide 5 of the ppt file.

## Data spreadsheets and files description

### Data taken with Rig 1 for positive currents while ramping current up.

 The central point (step # 27) was 
at X=-1.746cm, Y = -3.670cm and Z= -3.413cm 
relative to the magnet center ( need to clarify if it is mechanical center?)
The Z position for each point is given in the spreadsheets.

A set of measurements were taken at the central point (#27) for a series of positive currents. 
* The results are in file: HB-central field.xls
* The data extracted in file: HB-central-B-I.dat
* The field measurements need to be corrected for the non-linearity error of the hall probe.
* The TOSCA calculations for the probe location are given in the spreadsheet.


### Data taken with Rig 1 for positive currents while ramping current up.

* The hall probe could be moved along the central axis in fixed steps.
* They were 57 step points. Measurements were done at positive currents of 1200A, 2000A, 3000A, 3500A and 4000A
* The results are in file: HB-CenterLine-LakeshoreProbe.xls
* The data and TOSCA extracted in files: Central_zscan_files/central_zscan_{current}A_pos.dat
* The field measurements need to be corrected for the non-linearity error of the hall probe.
* The TOSCA calculations for the probe location are given in the spreadsheet.



### Data taken with Rig 1 for hysterisis curves for ramping up positive and negative currents.

Data was taken  at location step #27. Data was taken
for  negative  currents with the current ramped up
and the down. Then data was taken with negative current ramping up and the magnet quenched.
Then data was taking for the magnet with positive current and ramping the up and down.

  * The results are in file: hysteresis curve 2-26-16 Pos-Neg.xlsx
  * The data was extracted for the positive current ramp up: pos-up-current.dat
  * The data was extracted for the positive current ramp down: pos-down-current.dat
  * The data was extracted for the negative current ramp up: neg-up-current.dat
  * The data was extracted for the negative current ramp down: neg-down-current.dat

### Data taken with Rig 2 at points along Z for different central fields at 5 XY locations with negative polarity .

Data taken at 3899.35 +/- 0.15A, 3499.50 +/- 0.15A, 2999.65 +/- 0.15A and 1999.90 +/- 0.15A

  * The results are in file: Mapping Central Field 2000-3900a Neg.xls
  * Extracted data for each current and hole ( hole==XY location) in files: Central_zcan/central_zscan_{current}a_hole#.dat

### Root Code to analyze spreadsheet data
Root script plot_routines.C the following routines:
1. Get_two_col_data_from_file for reading in two column comma separated data into two vectors.
2. Set_probe which uses vectors read-in from the Lakeshore-probe.dat and makes plots. 
   * It fits the positive field data with 6th order polynomial (TF1 pos_field_corr) and 
   * the negative field data with a 8th order polynomial (TF1 neg_field_corr). 
3. Double_t field_corr(Double_t field ) which takes the field and determines field correction. 
   * For positive fields, evaluates pos_field_corr and multiplies by -1.
   * For negative fields, evaluates neg_field_corr and multiplies by +1.
   * This way new_field = field + field_corr(Double_t field ) regardless of the sign of the field.
4. Plot_hb_bi(vector<double> & ,vector<double> & ) plots vectors read in from HB-central-B-I.dat. The TOSCA data is hardcode. 
5. void run_code() which reads in data from probe and HB-central-B-I.dat , calls Set_probe and then calls Plot_hb_bi
6. void Plot_hb_pos_neg_diff() which reads in data from probe and calls Set_probe. 
   * Reads in data from the hystersis study (neg-down-current.dat,neg-up-current.dat,pos-up-current.dat,pos-down-current.dat) and makes plots.
7. void run_plot_zscan_pos(TString current)
   * Calls Set_probe 
   * Reads in data from positive current zscan with files Central_zscan_files/central_zscan_"+current+"A_pos.dat
   * Call plot_zscan_pos which corrects the data for the probe non-linearity and then plots the data and tosca.
   * Calculates the integral, max and effective field length for data and TOSCA. Prints them to screen.
   * Results saved to file Central_zscan_files/central_zscan_neg_results.dat
8. void run_plot_zscan(TString c, TString n);
   * Calls Set_probe
   * Reads in data from negative current zscan with files Central_zscan_files/central_zscan_+"c"+A_hole+"n".dat
   * Call plot_zscan which corrects the data for the probe non-linearity and then plots the data and tosca.
   * Calculates the integral, max and effective field length for data and TOSCA. Prints them to screen.
   * Results saved to file Central_zscan_files/central_zscan_neg_results.dat
9. void plot_zscan_results()
   * Reads in Central_zscan_files/central_zscan_neg_results.dat
   * Reads in Central_zscan_files/central_zscan_pos_results.dat
   * Reads in Central_zscan_files/tosca_HB-AB-cases.dat  
   * Makes plots comparing Max B/I, Int B/I and Eff length versus I.



 



