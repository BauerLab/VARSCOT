#!/usr/bin/python
import os
import Queue
import subprocess
import sys
import threading


INPUT_VCF = '../ram/input/input.vcf'
INPUT_BED = '../ram/input/targets.bed'
OUTPUT_TEMPLATE = '../output/{}.txt'
FASTA = '../ram/input/genome.fa'
INDEX_PREFIX = '../ram/input/index/prefix'
TEMP_SPACE = '../ram/tmp'
MISMATCHES = 4

THREADS = 48
OVERWRITE = True
RETRY_ON_ERROR = True

with open(INPUT_VCF, 'r') as vcf_file:
    for line in vcf_file:
        # Find the last line of the header
        if line.startswith('#CHROM'):
            headings = line.rstrip('\r\n').split('\t')
            samples = headings[9:][:96]
            break


def worker(queue_object):
    me = threading.current_thread().name
    while True:
        try:
            index, sample = queue_object.get(block=False)
        except Queue.Empty:
            sys.stdout.write(me + ' Queue empty\n')
            break
        sys.stdout.write(me + ' got ' + str((index, sample)) + '\n')
        temp_dir = os.path.join(TEMP_SPACE, str(index))
        worker_env = os.environ.copy()
        # Add temp dir to envs for bidir_index use
        worker_env['TMPDIR'] = temp_dir
        output_file = OUTPUT_TEMPLATE.format(sample)
        if os.path.exists(output_file) and not OVERWRITE:
            sys.stdout.write(me + ' Skipping\n')
            queue_object.task_done()
            continue
        try:
            output = subprocess.check_output(
                args=[
                    './VARSCOT',
                    '--vcf', INPUT_VCF,
                    '--bed', INPUT_BED,
                    '--output', output_file,
                    '--genome', FASTA,
                    '--index', INDEX_PREFIX,
                    '--sample', str(index),
                    '--mismatch', str(MISMATCHES),
                    # The whole process will run in one thread
                    '--threads', '1',
                    # pass temp dir to script for general use
                    '--temp-dir', temp_dir,
                ],
                env=worker_env,
            )
        except subprocess.CalledProcessError as e:
            sys.stdout.write(e.output)
            if RETRY_ON_ERROR:
                sys.stdout.write(me + ' Failed, adding back to queue\n')
                queue_object.put((index, sample))
            else:
                sys.stdout.write(me + ' Failed, skipping.\n')
            queue_object.task_done()
            continue
        sys.stdout.write(me + ' Completed.\n')
        queue_object.task_done()


q = Queue.Queue()

for index, sample in enumerate(samples):
    q.put((index, sample))

sys.stdout.write('Added ' + str(len(samples)) + ' samples\n')

threads = []
for i in range(THREADS):
    t = threading.Thread(target=worker, args=(q,))
    t.start()
    threads.append(t)

sys.stdout.write('Started ' + str(len(threads)) + ' threads\n')
# block until all tasks are done
q.join()

sys.stdout.write('Queue completed.\n')
# stop workers
for t in threads:
    t.join()

sys.stdout.write('All threads joined.\n')
