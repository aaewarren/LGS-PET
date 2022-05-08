#!/bin/bash

#Authored by Aaron Warren (aaron.warren@unimelb.edu.au)

#This code documents the image pre-processing steps relevant to the following paper:

#Frontoparietal 18FDG-PET hypo-metabolism in Lennox- Gastaut syndrome: further evidence for a shared epileptic network
#Tom Balfroid, Aaron E.L Warren, Linda J. Dalic, Alec Aeby, Salvatore U. Berlangieri, John S. Archer

#Analysis performed using the High Performance Computing (HPC) system ("Spartan") operated by Research Computing Services at The University of Melbourne:

#Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa
#https://dashboard.hpc.unimelb.edu.au 

#This bash script creates slurm submission scripts for each subject (*_sub.script). After per-subject scripts are created, they are ready to submit via sbatch.

#Usage: bash LGS_PET_PP.sh

#For each subject, the following steps are performed:

#1. Generate brain and white-matter masks from T1-weighted MRI using FreeSurfer's recon-all 
#2. Warp PET image to MNI symmetric template space by calculating and then combining T1-to-MNI and PET-to-T1 registrations
#3. If required, perform hemisphere flipping/deletion (remove hemispheres corresponding to side of epilepsy in TLE controls, and remove hemispheres with structural abnormalities in LGS patients, then flip kept hemispheres to same [right] side)
#4. Intensity normalise the resulting unilateral PET image using the mean signal within white matter (calculated from the same image side only)
#5. Apply spatial smoothing (FWHM=8mm) using FSL's SUSAN 

#After these pre-processing steps, the data are submitted to a between-group comparison (LGS vs TLE controls) using FSL's Permutation Analysis of Linear Models (PALM) software. See manuscript for details concerning this. 

#Load in a three-column .txt file called "SUBJECTS_GROUPS_HEMFLIP.txt" listing subject IDs (e.g., 001, 002, etc) in column 1, the group they belong to (LGS, L-TLE, R-TLE) in column 2, and whether or not hemisphere flipping
#is performed (0=no; 1=yes) in column 3

SUBJECTS=$(cat /data/PET/SUBJECTS_GROUPS_HEMFLIP.txt | awk '{print $1}'); #first column lists subject IDs

#Note the script assumes the PET and T1 MRI data are already in .nii.gz format and saved within per-subject sub-directories
#T1 should be called ${SUBJID}_MRI.nii.gz and be in ${anatdir}; PET should be called ${SUBJID}_PET.nii.gz and be in ${petdir}

#create slurm job for each subject
for SUBJID in $SUBJECTS; do 
	
	SCRIPTNAME=$SUBJID"_sub.script" 
	echo '#!/bin/bash' > $SCRIPTNAME
	echo '#SBATCH --time=36:0:0' >> $SCRIPTNAME
	echo '#SBATCH --ntasks=1' >> $SCRIPTNAME
	echo '#SBATCH --mem-per-cpu=32000' >> $SCRIPTNAME
	echo "#SBATCH --job-name=$SUBJID" >> $SCRIPTNAME
	echo >> $SCRIPTNAME 
		
	#load the software (+ build locations/setup scripts/dependencies, where necessary) we need
	echo module load gcc/8.3.0 >> $SCRIPTNAME
	echo module load openmpi/3.1.4 >> $SCRIPTNAME
	echo module load ants/2.3.2-python-3.7.4 >> $SCRIPTNAME #ANTs
	echo module load fsl/6.0.1-python-3.7.4 >> $SCRIPTNAME #FSL 
	echo module load mrtrix/3.0_rc2-python-2.7.16 >> $SCRIPTNAME #MRTrix 
	echo module load freesurfer/7.1.1-centos7_x86_64 >> $SCRIPTNAME	#FreeSurfer
	echo source /usr/local/easybuild/software/FSL/6.0.1-spartan_gcc-8.1.0-Python-3.7.3/fsl/etc/fslconf/fsl.sh >> $SCRIPTNAME
	
	echo SUBJECTS_DIR=/data/PET/FREESURFER_V7.1.1 >> $SCRIPTNAME #SUBJECTS_DIR where FreeSurfer recon-all output will be saved
	echo source /usr/local/easybuild-2019/easybuild/software/core/freesurfer/7.1.1-centos7_x86_64/SetUpFreeSurfer.sh >> $SCRIPTNAME
	echo MNI_DIR=/usr/local/easybuild-2019/easybuild/software/core/freesurfer/7.1.1-centos7_x86_64/mni >> $SCRIPTNAME
	
	#create/specify some directories to save output	(and where input MRI/PET files are)
	parentdir=/data/PET/${SUBJID};
	anatdir=${parentdir}/MRI
	petdir=${parentdir}/PET
	regdir=${parentdir}/REG2MNISYMM
	wmseg=${anatdir}/wmseg
	templatedir=/data/PET/MNI_152_6THGEN_SYMMETRIC #where the MNI template images and brain/hemisphere masks are. templates obtained from here: http://nist.mni.mcgill.ca/mni-icbm152-non-linear-6th-generation-symmetric-average-brain-stereotaxic-registration-model/

	echo "if [ ! -d ${wmseg} ]; then mkdir ${wmseg}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${regdir} ]; then mkdir ${regdir}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${SUBJECTS_DIR} ]; then mkdir ${SUBJECTS_DIR}; fi" >> $SCRIPTNAME
	
	#run freesurfer
	echo recon-all -i ${anatdir}/${SUBJID}_MRI.nii.gz -subjid $SUBJID -all >> $SCRIPTNAME 

	#convert freesurfer version of T1 to nii.gz format using mrtrix
	echo mrconvert ${SUBJECTS_DIR}/${SUBJID}/mri/orig/001.mgz ${anatdir}/origanat.nii.gz -stride +1,+2,+3 -force >> $SCRIPTNAME
	
	#create head mask (not brain mask) of anat file using bet
	echo bet $anatdir/origanat.nii.gz $anatdir/bet -A >> $SCRIPTNAME
	
	#delete bet files we don't need
	echo rm -f $anatdir/bet*mesh* $anatdir/bet.nii.gz $anatdir/bet*skull* >> $SCRIPTNAME
	
	#rename the head mask
	echo mv $anatdir/bet_outskin_mask.nii.gz $anatdir/anat_headmask.nii.gz >> $SCRIPTNAME
	
	#generate masked head
	echo fslmaths $anatdir/origanat.nii.gz -mul $anatdir/anat_headmask.nii.gz $anatdir/anat_head.nii.gz >> $SCRIPTNAME
	
	#bias correct the head-masked T1 using ANTs
	echo N4BiasFieldCorrection -d 3 -i ${anatdir}/anat_head.nii.gz -x ${anatdir}/anat_headmask.nii.gz -o ${anatdir}/anat_head_biascorr.nii.gz >> $SCRIPTNAME
	
	#create a brain mask (at resolution of input T1) using the output from freesurfer (aparc+aseg parcellation), then fill any holes, and lastly resample everything to 09mm iso res 
	echo mri_label2vol --seg ${SUBJECTS_DIR}/${SUBJID}/mri/aparc+aseg.mgz --temp ${SUBJECTS_DIR}/${SUBJID}/mri/rawavg.mgz --o ${anatdir}/aparc+aseg2rawavg.mgz --regheader ${SUBJECTS_DIR}/${SUBJID}/mri/aparc+aseg.mgz >> $SCRIPTNAME
	echo mrconvert ${anatdir}/aparc+aseg2rawavg.mgz ${anatdir}/aparc+aseg.nii.gz -stride +1,+2,+3 -force >> $SCRIPTNAME
	echo fslmaths ${anatdir}/aparc+aseg.nii.gz -bin -kernel sphere 3 -dilM -ero -dilM -ero -dilM -ero -binv ${anatdir}/anat_mask_tmp1.nii.gz >> $SCRIPTNAME
	echo connectedcomp ${anatdir}/anat_mask_tmp1.nii.gz ${anatdir}/anat_mask_tmp2.nii.gz >> $SCRIPTNAME
	echo fslmaths ${anatdir}/anat_mask_tmp2.nii.gz -uthr 1 -binv ${anatdir}/anat_brainmask.nii.gz >> $SCRIPTNAME
	echo rm ${anatdir}/anat_mask_tmp?.nii.gz ${anatdir}/aparc+aseg2rawavg.mgz >> $SCRIPTNAME
	echo fslmaths ${anatdir}/anat_head_biascorr.nii.gz -mul ${anatdir}/anat_brainmask.nii.gz ${anatdir}/anat_brain_biascorr.nii.gz >> $SCRIPTNAME
	echo flirt -in ${anatdir}/anat_head_biascorr.nii.gz -ref ${anatdir}/anat_head_biascorr.nii.gz -out ${anatdir}/anat_head_biascorr_09mm.nii.gz -nosearch -applyisoxfm 0.9 >> $SCRIPTNAME
	echo flirt -in ${anatdir}/anat_brain_biascorr.nii.gz -ref ${anatdir}/anat_brain_biascorr.nii.gz -out ${anatdir}/anat_brain_biascorr_09mm.nii.gz -nosearch -applyisoxfm 0.9 >> $SCRIPTNAME
	echo fslmaths ${anatdir}/anat_head_biascorr_09mm.nii.gz -bin ${anatdir}/anat_headmask_09mm.nii.gz >> $SCRIPTNAME
	echo fslmaths ${anatdir}/anat_brain_biascorr_09mm.nii.gz -bin ${anatdir}/anat_brainmask_09mm.nii.gz >> $SCRIPTNAME
	echo flirt -in ${anatdir}/aparc+aseg.nii.gz -ref ${anatdir}/aparc+aseg.nii.gz -out ${anatdir}/aparc+aseg_09mm.nii.gz -nosearch -applyisoxfm 0.9 -interp nearestneighbour >> $SCRIPTNAME
	
	#warp to symmetric MNI 6th gen space (1mm iso) using ANTs, with both whole-head and brain-extracted images driving the registration
	echo antsRegistration --verbose 1 --dimensionality 3 --float 0 --output [${regdir}/anat_to_mni_symm_1mm,${regdir}/anat_to_mni_symm_1mm.nii.gz] \
		--interpolation Linear --use-histogram-matching 1 --winsorize-image-intensities [0.005,0.995] --write-composite-transform 1 \
		--initial-moving-transform [${templatedir}/MNI152_T1_1mm_symmetric.nii.gz,${anatdir}/anat_head_biascorr_09mm.nii.gz,1] \
		--transform Rigid[0.15] --metric MI[${templatedir}/MNI152_T1_1mm_symmetric.nii.gz,${anatdir}/anat_head_biascorr_09mm.nii.gz,1,32,Regular,0.25] \
		--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.15] \
		--metric MI[${templatedir}/MNI152_T1_1mm_symmetric.nii.gz,${anatdir}/anat_head_biascorr_09mm.nii.gz,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] \
		--shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[0.1,3,0] --metric CC[${templatedir}/MNI152_T1_1mm_symmetric.nii.gz,${anatdir}/anat_head_biascorr_09mm.nii.gz,0.3,4] \
		--metric CC[${templatedir}/MNI152_T1_1mm_brain_symmetric.nii.gz,${anatdir}/anat_brain_biascorr_09mm.nii.gz,0.7,4] --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 \
		--smoothing-sigmas 3x2x1x0vox >> $SCRIPTNAME

	#register pet to t1 using ANTs
	echo antsRegistration --dimensionality 3 --float 0 --output [${petdir}/PET_to_anat_09mm,${petdir}/PET_to_anat_09mm.nii.gz] --interpolation Linear \
		--use-histogram-matching 0 --initial-moving-transform [${anatdir}/anat_brain_biascorr_09mm.nii.gz,${petdir}/${SUBJID}_PET.nii.gz,1] --transform Rigid[0.1] \
		--metric MI[${anatdir}/anat_brain_biascorr_09mm.nii.gz,${petdir}/${SUBJID}_PET.nii.gz,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 \
		--smoothing-sigmas 3x2x1x0vox >> $SCRIPTNAME
	
	#warp pet to symmetric mni space (1mm iso) using non-linear and linear transforms calculated above
	echo antsApplyTransforms --dimensionality 3 --float 0 --input ${petdir}/${SUBJID}_PET.nii.gz --reference-image ${templatedir}/MNI152_T1_1mm_symmetric.nii.gz \
		--output ${regdir}/PET_to_mni_symm_1mm.nii.gz --interpolation Linear --transform ${regdir}/anat_to_mni_symm_1mmComposite.h5 --transform ${petdir}/PET_to_anat_09mm0GenericAffine.mat >> $SCRIPTNAME
		
	#generate white-matter mask using recon-all output, by firstly extracting desired white-matter segments from aparc+aseg. 
	echo "for roi in Left-Cerebral-White-Matter Right-Cerebral-White-Matter CC_Posterior CC_Mid_Posterior CC_Central CC_Mid_Anterior CC_Anterior; do" >> $SCRIPTNAME

		#get code from the freesurfer lookup table
		echo 'roi_code=`cat /usr/local/easybuild-2019/easybuild/software/core/freesurfer/7.1.1-centos7_x86_64/FreeSurferColorLUT.txt | grep -w ${roi}" " | awk '"'"'{print $1}'"'"'`' >> $SCRIPTNAME #note the complicated single/double quotes to prevent bash eval 
		echo 'thr_lower=`echo "scale=1; ${roi_code} - 0.1" | bc `' >> $SCRIPTNAME
		echo 'thr_upper=`echo "scale=1; ${roi_code} + 0.1" | bc `' >> $SCRIPTNAME
		
		#extract the roi from aparc+aseg segmentation 
		echo 'fslmaths '${anatdir}'/aparc+aseg_09mm.nii.gz -thr ${thr_lower} -uthr ${thr_upper} -bin '${wmseg}'/roi_t1_${roi}.nii.gz' >> $SCRIPTNAME
		
	echo "done" >> $SCRIPTNAME

	#add all the white matter segments together to make one big mask 
	echo fslmaths ${wmseg}/roi_t1_Left-Cerebral-White-Matter.nii.gz -add ${wmseg}/roi_t1_Right-Cerebral-White-Matter.nii.gz \
		-add ${wmseg}/roi_t1_CC_Posterior.nii.gz -add ${wmseg}/roi_t1_CC_Mid_Posterior.nii.gz -add ${wmseg}/roi_t1_CC_Central.nii.gz \
		-add ${wmseg}/roi_t1_CC_Mid_Anterior.nii.gz -add ${wmseg}/roi_t1_CC_Anterior.nii.gz -bin ${wmseg}/roi_t1_Freesurfer-Cerebral-WM.nii.gz >> $SCRIPTNAME
	
	#warp the white matter mask into 1mm symmetric mni space, using trilinear interpolation 
	echo antsApplyTransforms --dimensionality 3 --float 0 --input ${wmseg}/roi_t1_Freesurfer-Cerebral-WM.nii.gz --reference-image ${templatedir}/MNI152_T1_1mm_symmetric.nii.gz \
		--output ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_prob.nii.gz --interpolation Linear --transform ${regdir}/anat_to_mni_symm_1mmComposite.h5 >> $SCRIPTNAME
	
	#threshold at 0.99
	echo fslmaths ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_prob.nii.gz -thr 0.99 -bin ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin.nii.gz >> $SCRIPTNAME
	
	#determine whether hemisphere flipping is required
	#note: flipping is only performed for patients with RIGHT sided abnormalities - i.e., R-TLE, or LGS patients with right hem structural abnormalities. the brain is flipped so that the left hemisphere appears on the right, and then the "was right but now appears on left hem" is deleted. 	
	#in contrast, for patients with left-sided abnormalities, no flipping is performed (i.e., the right hemisphere is kept in place and the left hemisphere is just deleted)

	flip=$(cat /data/PET/SUBJECTS_GROUPS_HEMFLIP.txt | grep ${SUBJID} | awk '{print $3}'); #1 = yes do flipping; 0 = don't do any flipping
	
	if [ "$flip" = "1" ]; then #yes do flipping
		
		#flip about x axis 
		echo fslswapdim ${regdir}/PET_to_mni_symm_1mm.nii.gz -x y z ${regdir}/PET_to_mni_symm_1mm_LEFTISNOWRIGHT.nii.gz >> $SCRIPTNAME

		#remove the new "left" (originally right) hemisphere, preserving the new "right" (originally left) hemisphere
		echo fslmaths ${regdir}/PET_to_mni_symm_1mm_LEFTISNOWRIGHT.nii.gz -mul ${templatedir}/MNI152_T1_1mm_brain_symmetric_rh_brainmask.nii.gz ${regdir}/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked.nii.gz >> $SCRIPTNAME
		
		#also flip the wm segmentation in 1mm symm MNI space 	
		echo fslswapdim ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin.nii.gz -x y z ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_LEFTISNOWRIGHT.nii.gz >> $SCRIPTNAME
		
		#and mask with rh mask 
		echo fslmaths ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_LEFTISNOWRIGHT.nii.gz -mul ${templatedir}/MNI152_T1_1mm_brain_symmetric_rh_brainmask.nii.gz ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_LEFTISNOWRIGHT_rhmasked.nii.gz >> $SCRIPTNAME
	
		#erode white matter mask by zeroing non-zero voxels when zero voxels are found within a 5mm radius sphere centred on each voxel
		echo fslmaths ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_LEFTISNOWRIGHT_rhmasked.nii.gz -kernel sphere 5 -ero ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_LEFTISNOWRIGHT_rhmasked_ero5mmradsphere.nii.gz >> $SCRIPTNAME
				
		#WHITE MATTER INTENSTIY NORMALISATION#

		#calculate mean PET signal within "right" (originally left) hemsiphere eroded white matter mask created above
		echo 'wm_mean=$(fslstats '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked.nii.gz -k '${wmseg}'/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_LEFTISNOWRIGHT_rhmasked_ero5mmradsphere.nii.gz -m)' >> $SCRIPTNAME

		#now do intensity normalisation: divide by right wm mean
		echo 'fslmaths '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked.nii.gz -div $wm_mean '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked_wmnorm.nii.gz' >> $SCRIPTNAME

		#smooth (only right hemisphere) at 8mm and re-mask, using SUSAN, and resample to 2mm for analysis
		fwhm=8
		echo 'thresh=$(fslstats '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked_wmnorm.nii.gz -k ${templatedir}/MNI152_T1_1mm_brain_symmetric_rh_brainmask.nii.gz -p 50)' >> $SCRIPTNAME #extracts 50th-percentile of values within pet image			
		echo 'thresh75="$(echo "$thresh*0.75" | bc)"' >> $SCRIPTNAME #find 75% of the 50th-percentile calculate in previous line: this is default threshold recommended for susan	
		echo 'susan '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked_wmnorm.nii.gz $thresh75 '${fwhm}' 3 1 1 '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked_wmnorm.nii.gz $thresh75 '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz' >> $SCRIPTNAME
		echo 'flirt -in '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -ref '${regdir}'/PET_to_mni_symm_1mm_LEFTISNOWRIGHT_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -out '${regdir}'/PET_to_mni_symm_2mm_LEFTISNOWRIGHT_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -nosearch -applyisoxfm 2.0' >> $SCRIPTNAME	
		
		#re-mask 			
		echo 'fslmaths '${regdir}'/PET_to_mni_symm_2mm_LEFTISNOWRIGHT_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -mul ${templatedir}/MNI152_T1_2mm_brain_symmetric_rh_brainmask.nii.gz '${regdir}'/PET_to_mni_symm_2mm_LEFTISNOWRIGHT_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz' >> $SCRIPTNAME
		
	elif [ "$flip" = "0" ]; then #no flipping
		
		#remove the left hemisphere, preserving the right
		echo fslmaths ${regdir}/PET_to_mni_symm_1mm.nii.gz -mul ${templatedir}/MNI152_T1_1mm_brain_symmetric_rh_brainmask.nii.gz ${regdir}/PET_to_mni_symm_1mm_rhmasked.nii.gz >> $SCRIPTNAME
		
		#also rh mask the wm segmentation in 1mm symm MNI space  
		echo fslmaths ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin.nii.gz -mul ${templatedir}/MNI152_T1_1mm_brain_symmetric_rh_brainmask.nii.gz ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_rhmasked.nii.gz >> $SCRIPTNAME
	
		#erode white matter mask by zeroing non-zero voxels when zero voxels are found within a 5mm radius sphere centred on each voxel
		echo fslmaths ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_rhmasked.nii.gz -kernel sphere 5 -ero ${wmseg}/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_rhmasked_ero5mmradsphere.nii.gz >> $SCRIPTNAME
				
		#WHITE MATTER INTENSTIY NORMALISATION#

		#calculate mean PET signal within "right" (originally left) hemsiphere eroded white matter mask created above
		echo 'wm_mean=$(fslstats '${regdir}'/PET_to_mni_symm_1mm_rhmasked.nii.gz -k '${wmseg}'/roi_mni_symm_1mm_Freesurfer-Cerebral-WM_thresh99_bin_rhmasked_ero5mmradsphere.nii.gz -m)' >> $SCRIPTNAME

		#now do intensity normalisation: divide by right wm mean
		echo 'fslmaths '${regdir}'/PET_to_mni_symm_1mm_rhmasked.nii.gz -div $wm_mean '${regdir}'/PET_to_mni_symm_1mm_rhmasked_wmnorm.nii.gz' >> $SCRIPTNAME

		#smooth (only right hemisphere) at 8mm and re-mask, using SUSAN, and resample to 2mm for analysis
		fwhm=8
		echo 'thresh=$(fslstats '${regdir}'/PET_to_mni_symm_1mm_rhmasked_wmnorm.nii.gz -k ${templatedir}/MNI152_T1_1mm_brain_symmetric_rh_brainmask.nii.gz -p 50)' >> $SCRIPTNAME #extracts 50th-percentile of values within pet image		
		echo 'thresh75="$(echo "$thresh*0.75" | bc)"' >> $SCRIPTNAME #find 75% of the 50th-percentile calculate in previous line: this is default threshold recommended for susan			
		echo 'susan '${regdir}'/PET_to_mni_symm_1mm_rhmasked_wmnorm.nii.gz $thresh75 '${fwhm}' 3 1 1 '${regdir}'/PET_to_mni_symm_1mm_rhmasked_wmnorm.nii.gz $thresh75 '${regdir}'/PET_to_mni_symm_1mm_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz' >> $SCRIPTNAME
		echo 'flirt -in '${regdir}'/PET_to_mni_symm_1mm_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -ref '${regdir}'/PET_to_mni_symm_1mm_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -out '${regdir}'/PET_to_mni_symm_2mm_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -nosearch -applyisoxfm 2.0' >> $SCRIPTNAME
			
		#re-mask 			
		echo 'fslmaths '${regdir}'/PET_to_mni_symm_2mm_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz -mul ${templatedir}/MNI152_T1_2mm_brain_symmetric_rh_brainmask.nii.gz '${regdir}'/PET_to_mni_symm_2mm_rhmasked_wmnorm_susan'${fwhm}'mm.nii.gz' >> $SCRIPTNAME
			
	fi
	
done