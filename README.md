![Intro](/Commersons-Peales-hydrophone.jpg)

# Introduction
This project is part of a sound engineering thesis, treating a specific need brought up by Fundación Cethus (NGO - Argentina).

The implementation of acoustic methods used to study cetaceans is increasing due to their economic advantages in comparison to visual monitoring methods. These acoustic methods present challenges both during recording stage and post processing. One particular case is the classification of similar acoustic signals emitted by sympatric dolphin species. The high directivity of echolocation signals as well as the side to side movements of the head of animals while scanning makes the recording process and characterization of such signals a difficult task. Classification of NBHF clicks from sympatric species is a necessity in order to study them through acoustic methods. 

A MATLAB code was developed in order to detect Peale´s and Commerson’s dolphins’ clicks, characterize them, and classify them according to the species which produced them. A semi-automatic detector was designed to find peaks by setting a 24 dB SNR threshold. Both on and off-axis dolphins’ clicks were extracted, separating them from unwanted noise from a set of audio files. Peak frequency, centroid frequency, -3 dB and -10 dB bandwidths, and -10 dB duration were obtained for each detected click. A ‘Support Vector Machines’ algorithm was implemented for the classification process, using the free software package LIBSVM. Several evaluations were carried out while searching for the optimum performance of a class prediction model capable of classifying clicks according to emitting dolphin species. 

It is concluded that clicks produced by Peale´s and Commerson’s dolphins in Santa Cruz province in Argentina can be separated with an estimate precision of 90% considering a series of future considerations to enhance and further develop the current work.

# Methodology
## Data collection
Data was collected in two different ways:
- By using a towed hydrophone array from a large ship
- By using a dip hydrophone from a small boat

Acoustic recordings were sammpled at 500 kHz and 16-bit resolution into 5 minute files.
Three sets of data were gathered:

- Only Commerson's dolphins present in the area
- Only Peale's dolphins present in the area
- Both species producing clicks simultaneously

Presence of specific species were visually confirmed by a dedicated observer during recordings.

## Data processing
### Detection
Long Term Spectral Averages (LTSAs) were visually inspected for a first discrimination of files likely to include NBHF clicks from dolphins. The selected files were looped over and analysed by the main script **DetectorClicks.m** which includes calls to **detector_umbral.m** for click detection witihin each file and **acoustic_params.m** for acoustic signal description of each stored click. The output of the script is a structure with the following format:
<pre>

├───File (1)
|  ├─── Det (1)
|  |     ├─── filename (where the click can be found)
|  |     ├─── clickn (detection number within file)
|  |     ├─── itime (initial click time referenced to beginning of file - datestr format)
|  |     |─── spectrum (magnitude spectrum of the click)
|  |     ├─── signal (time series amplitud values for the detection window of the click)
|  |     ├─── acoustic_params (time and spectral domain signal features)
|  |
|  ├─── Det (2)
|  :
|
├─── File (2)  
:

</pre>
This structure allows to locate and go back to any specific click or click train to view in detail possible causes of outliers, malfunctions in the detection function, or to review in detail its waveform and spectrum. 

### Classification
Information used to train and validate the **support vector machine** model is stored in a matrix called "SVMmat", which is then passed to the **ClasificadorClicks.m** script. LibSVM library was used for the SVM implementation. Validation was carried out using *5-fold Cross Validation* and wightening used to correct for unbalanced dataset sizes.

After the model was built, it was tested using the recordings where both species are suposed to be overlapping their click trains.

# Results
For visualization purposes a small portion of one of the files where both species were present was plotted, showing red and blue dots corresponding to times of detected NBHF clicks and each color displaying the output of the classifier for that specific click.

![Results](/Classifier-output.jpg)

# Publications

Results were presented in:

1.	**Marino, A.**, Reyes Reyes, M.V., & Melcón, M.L. Acoustic detection and classification of NBHF clicks from free-ranging Commerson´s and Peale´s dolphins. Accepted for presentation at the 3rd Listening for Aquatic Mammals in Latin America Workshop. Lima, Perú, 4th November 2018. (Oral presentation).

2.	**Marino, A.**, Reyes Reyes, M.V., & Melcón, M.L. Detección acústica y clasificación de clicks NBHF de delfines australes y toninas overas en su hábitat natural. Accepted for presentation at the XVIII Conferences on Specialists on Aquatic Mammals from South America and 12th SOLAMAC Congress. Lima, Perú, 5th – 8th November 2018. (Poster).
