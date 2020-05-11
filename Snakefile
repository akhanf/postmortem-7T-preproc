from os.path import join
from glob import glob

configfile: "config.yml"


# bids_dir  set in json file.
# can also override at command-line with e.g.:  --config bids_dir='path/to/dir'  or --configfile ...
bids_dir = config['bids_dir']
subjects = config['subjects']


rule all:
    input:
        mni_affine = expand('sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct_to-MNI152NLin2009cAsym_affine.nii.gz',subject=subjects),
        mni_rigid = expand('sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct_to-MNI152NLin2009cAsym_rigid.nii.gz',subject=subjects)


rule align_run_flirt:
    input:
        fixed = lambda wildcards: join(bids_dir,'sub-{subject}/anat/sub-{subject}_acq-SPACE_run-0{run}_T2w.nii.gz'.format(subject=wildcards.subject,run=config['use_runs'][wildcards.subject][0])),
        moving = join(bids_dir,'sub-{subject}/anat/sub-{subject}_acq-SPACE_run-0{run}_T2w.nii.gz')
    output:
        warped = 'sub-{subject}/sub-{subject}_run-0{run}_T2w_regToRef_warped.nii.gz',
        xfm = 'sub-{subject}/sub-{subject}_run-0{run}_T2w_regToRef_xfm.mat' 
    envmodules: 'fsl'
    log: 'logs/align_run_flirt/sub-{subject}_run-0{run}.log'
    shell:
        'flirt -in {input.moving} -ref {input.fixed} -out {output.warped} -omat {output.xfm} -dof 6 -coarsesearch 30 -finesearch 15 &> {log}'
 
rule merge_runs:
    input: 
        aligned = lambda wildcards: expand('sub-{subject}/sub-{subject}_run-0{run}_T2w_regToRef_warped.nii.gz',subject=wildcards.subject,run=config['use_runs'][wildcards.subject][1:]),
        fixed = lambda wildcards: join(bids_dir,'sub-{subject}/anat/sub-{subject}_acq-SPACE_run-0{run}_T2w.nii.gz'.format(subject=wildcards.subject,run=config['use_runs'][wildcards.subject][0])),
    output:
        aligned_4d = 'sub-{subject}/sub-{subject}_T2w_aligned_4d.nii.gz'
    envmodules: 'fsl'
    log: 'logs/merge_runs/sub-{subject}.log'
    shell:
        'fslmerge -t {output} {input} &> {log}'

rule avg_runs:
    input:
        aligned_4d = 'sub-{subject}/sub-{subject}_T2w_aligned_4d.nii.gz'
    output:
        avg = 'sub-{subject}/sub-{subject}_T2w_aligned_avg.nii.gz' 
    envmodules: 'fsl'
    log: 'logs/avg_runs/sub-{subject}.log'
    shell:
        'fslmaths {input} -Tmean {output}  &> {log}'

#run n4 multiple times, based on num_n4_reps['subject'] 
rule nu_correct:
    input: 
        avg = 'sub-{subject}/sub-{subject}_T2w_aligned_avg.nii.gz' 
    params:
        n4_reps = lambda wildcards: config['num_n4_reps'][wildcards.subject]
    output: 
        n4_dir = directory('sub-{subject}/sub-{subject}_T2w_aligned_avg_n4iters'),
        corr = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct.nii.gz',
        biasfield = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4biasfield.nii.gz'
    envmodules: 'ants'
    log: 'logs/nu_correct/sub-{subject}.log'
    shell:
        'mkdir -p {output.n4_dir} && ln -srv {input} {output.n4_dir}/n4_corr_0.nii.gz  &> {log} && ' 
        'for i in `seq 0 {params.n4_reps}`; do '
        '  N4BiasFieldCorrection -i {output.n4_dir}/n4_corr_$i.nii.gz -o [{output.n4_dir}/n4_corr_$((i+1)).nii.gz,{output.n4_dir}/n4_biasfield_$((i+1)).nii.gz] &>> {log}; '
        'done &&'
        'cp -v {output.n4_dir}/n4_corr_{params.n4_reps}.nii.gz {output.corr} &&'
        'cp -v {output.n4_dir}/n4_biasfield_{params.n4_reps}.nii.gz {output.biasfield}'

#register to mni (T2w brain) template
rule mni_reg_flirt_affine:
    input:
        moving = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct.nii.gz',
        fixed = 'template/tpl-MNI152NLin2009cAsym_res-01_desc-brain_T2w.nii.gz',
        tform =  lambda wildcards: 'sub-{subject}/sub-{subject}_reorient_tform.mat'.format(subject=wildcards.subject) if config['enable_custom_flirt'][wildcards.subject] else ''
    params:
        init = lambda wildcards, input: '-init {tform}'.format(tform=input.tform) if config['enable_custom_flirt'][wildcards.subject] else ''
    output:
        warped = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct_to-MNI152NLin2009cAsym_affine.nii.gz',
        xfm = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct_to-MNI152NLin2009cAsym_affine_xfm.mat'
    envmodules: 'fsl'
    log: 'logs/mni_reg_flirt/sub-{subject}.log'
    shell:
        'flirt -in {input.moving} -ref {input.fixed} -out {output.warped} {params.init} -omat {output.xfm} -dof 12 -coarsesearch 30 -finesearch 15 &> {log}'
 
#register to mni (T2w brain) template
rule mni_reg_flirt_rigid:
    input:
        moving = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct.nii.gz',
        fixed = 'template/tpl-MNI152NLin2009cAsym_res-01_desc-brain_T2w.nii.gz',
        tform =  lambda wildcards: 'sub-{subject}/sub-{subject}_reorient_tform.mat'.format(subject=wildcards.subject) if config['enable_custom_flirt'][wildcards.subject] else ''
    params:
        init = lambda wildcards, input: '-init {tform}'.format(tform=input.tform) if config['enable_custom_flirt'][wildcards.subject] else ''
    output:
        warped = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct_to-MNI152NLin2009cAsym_rigid.nii.gz',
        xfm = 'sub-{subject}/sub-{subject}_T2w_aligned_avg_n4correct_to-MNI152NLin2009cAsym_rigid_xfm.mat'
    envmodules: 'fsl'
    log: 'logs/mni_reg_flirt/sub-{subject}.log'
    shell:
        'flirt -in {input.moving} -ref {input.fixed} -out {output.warped} {params.init} -omat {output.xfm} -dof 12 -coarsesearch 30 -finesearch 15 &> {log}'
