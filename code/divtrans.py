import pysam
import logging 

logging.basicConfig(filename='divtrans.log',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

class Oread(object):
    """ An oriented read objet
    """
    def __init__(self, ali_read):
        """Extract only monolitic block positions (no splicing considered) and
        read orientation from a pysam alignment"""

        blocks = ali_read.get_blocks()
        self.chro = ali_read.reference_name
        self.start = blocks[0][0]
        self.end = blocks[-1][-1]

        if (ali_read.is_read1 and ali_read.is_reverse) or (ali_read.is_read2 and not ali_read.is_reverse):
            self.orientation = '+'
            
        elif (ali_read.is_read1 and not ali_read.is_reverse) or (ali_read.is_read2 and ali_read.is_reverse):
            self.orientation = '-'
        
        else:
            raise(IOError('Check loop !'))

    def __repr__(self):
        return "<Oread object: {chro}:{start}-{end} ({ori})>".format(chro=self.chro,
                                                                   start=self.start,
                                                                   end=self.end,
                                                                   ori=self.orientation)

def get_oriented_reads(bam):
    """To investigate the orientation of the transcription accross an alignment.
    Yield each read positions accross an alignment with its orientation
    """

    count = 0

    with pysam.AlignmentFile(bam) as bam:
        for read in bam:
        
            count += 1
            if count % 10000000 == 0:
                logging.debug("{} reads treated".format(count))
            
            yield Oread(read)

            
def identify_divergent_transcription(oread_gen, distance=100):
    """From a generator from get_transcription_strands, identify sites of
    divergent transcription.

    Distance specifies the maximum distance between divergent reads to
    consider a divergent transcription event
    """
    
    upstream_oread = next(oread_gen)

    for downstream_oread in oread_gen:

        if upstream_oread.orientation != '-':
            pass

        elif (downstream_oread.start - upstream_oread.end <= distance) and \
             (downstream_oread.orientation == '+'):

            if downstream_oread.chro == upstream_oread.chro:

                logging.debug("DIV found for reads {0} and {1}".format(upstream_oread,
                                                                       downstream_oread))

                yield (upstream_oread.chro,
                       min(upstream_oread.end, downstream_oread.start),
                       max(upstream_oread.end, downstream_oread.start))
            
        upstream_oread = downstream_oread


