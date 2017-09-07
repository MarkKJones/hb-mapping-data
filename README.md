#HB mapping data


A variety of mapping data was taken for HB by Eric Sun, Steve Lassiter 
and Mike Fowler from roughly Dec 2015-Mar 2016.

## Hall Probe info

A Lakeshore probe Lakeshore-XHGT-9060-55937-Q1 was used for the measurements.
The probe correction table is in Lakeshore-XHGT-9060-55937-Q1.pdf .
The data from the table has been put in file: Lakeshore-probe.dat
This linearity error is different for positive and negative field directions. 
The general Lakeshore manual is file: Lakeshore_general_manual.pdf
In section 4.5 in the Lakeshore manual it states the the error in a magnitude error.
In the data table, the Error = abs(Measured Field)-abs(True Field). 
So for positive fields True Field = Measured Field - Error. 
For negative fields True Field = Measured Field + Error.

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

The hall probe could be moved along the central axis in fixed steps.
They were 57 step points. The central point (step # 27) was 
at X=-1.746cm, Y = -3.670cm and Z= -3.413cm 
relative to the magnet center ( need to clarify if it is mechanical center?)
The Z position for each point is given in the spreadsheets.

A set of measurements were taken at the central point (#27) for a series of positive currents. 
The results are in file: HB-central field.xls
The data extracted in file: HB-central-B-I.dat
The field measurements need to be corrected for the non-linearity error of the hall probe.
The TOSCA calculations for the probe location are given in the spreadsheet.

### Data taken with Rig 1 for hysterisis curves for ramping up positive and negative currents.

Data was taken  at location step #27. Data was taken
for  negative  currents with the current ramped up
and the down. Then data was taken with negative current ramping up and the magnet quenched.
Then data was taking for the magnet with positive current and ramping the up and down.

The results are in file: hysteresis curve 2-26-16 Pos-Neg.xlsx
The data was extracted for the positive current ramp up: pos-up-current.dat
The data was extracted for the positive current ramp down: pos-down-current.dat
The data was extracted for the negative current ramp up: neg-up-current.dat
The data was extracted for the negative current ramp down: neg-down-current.dat





