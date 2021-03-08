import os
import subprocess
import tempfile
import logging

log = logging.getLogger('2passtools')

MINIMAP2 = os.environ.get('MINIMAP2_PATH', 'minimap2')
SAMTOOLS = os.environ.get('SAMTOOLS_PATH', 'samtools')

for program, prg_name in zip([MINIMAP2, SAMTOOLS], ['minimap2', 'samtools']):
    try:
        subprocess.check_call([program, '--help'])
    except FileNotFoundError:
        raise OSError(f'{prg_name} not found at location "{program}"')


def subprocess_command(cmd, stdout_fn):
    with open(stdout_fn, 'w') as s:
        proc = subprocess.Popen(
            cmd,
            stdout=s,
            stderr=subprocess.PIPE
        )
        _, stderr = proc.communicate()
        if proc.returncode:
            raise subprocess.CalledProcessError(stderr.decode())
        else:
            return stderr.decode().strip()


def map_with_minimap2(fastq_fn, reference_fn, output_fn=None, threads=1,
                      k=14, w=5, splice_flank=True, noncanon_pen=9,
                      canon_splice_strand='f',
                      junc_bed=None, junc_bonus=12,
                      max_intron_size=100000,
                      end_seed_pen=12):
    if not os.path.exists(fastq_fn):
        raise OSError('fastq_fn not found')
    elif not os.path.exists(reference_fn):
        raise OSError('reference_fn not found')
    splice_flank = 'yes' if splice_flank else 'no'
    s_handle, sam_fn = tempfile.mkstemp(suffix='.sam')
    b_handle, bam_fn = tempfile.mkstemp(suffix='.bam')

    # run minimap
    minimap2_cmd = [
        MINIMAP2, f'-t{threads}', '-a', '-x', 'splice',
        '-L', '--cs=long', f'-G{max_intron_size}',
        f'-C{noncanon_pen}', f'--splice-flank={splice_flank}',
        f'-u{canon_splice_strand}', f'--end-seed-pen={end_seed_pen}',
    ]
    if junc_bed is not None:
        minimap2_cmd += ['--junc-bed', junc_bed, f'--junc-bonus={junc_bonus}']
    minimap2_cmd += [reference_fn, fastq_fn]

    log.info('Running minimap2')
    log.debug('cmd: ' + ' '.join(minimap2_cmd))
    minimap2_stderr = subprocess_command(minimap2_cmd, sam_fn)
    log.info('Mapping complete')
    log.debug(f'minimap2 stderr:\n{minimap2_stderr}')

    # run samtools view
    samtools_view_cmd = [SAMTOOLS, 'view', '-bS', sam_fn]
    samtools_view_stderr = subprocess_command(samtools_view_cmd, bam_fn)
    

    # clean up minimap2 output
    os.close(s_handle)
    os.remove(sam_fn)

    # run samtools sort
    samtools_sort_cmd = [SAMTOOLS, 'sort', '-@', str(threads), '-o', '-', bam_fn]
    samtools_sort_stderr = subprocess_command(samtools_sort_cmd, output_fn)
    log.debug(f'samtools sort stderr:\n{samtools_sort_stderr}')

    # clean up samtools view output
    os.close(b_handle)
    os.remove(bam_fn)

    # run samtools index
    samtools_index_cmd = [SAMTOOLS, 'index', output_fn]
    subprocess.check_call(samtools_index_cmd)