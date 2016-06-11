from multiprocessing import Pool
import sys

class ProducerConsumer(object):
    def __init__(self, args, threads, producer, consumer):
        try:# in case we forget to set a debug flag as an option. set it up here for use..
          self.debug = args.debug
        except:
          print >> sys.stderr, "args is missing debug flag.  Setting it to False for ProducerConsumer"
          self.debug = False
        # we subtract 1 thread, to account for the consumer thread. aka the main thread
        if threads > 1:
            self.threads= threads - 1
        else:
            self.threads = 1
        self.producer = producer
        self.consumer = consumer

    def run(self, passthrough, producerargs):
        if not self.debug:
            pool = Pool(processes = self.threads)
            ret = pool.imap_unordered(self.producer, producerargs)
            self.consumer(passthrough, ret)
            pool.close()
        else:
            self.consumer(passthrough, ( self.debugwrap(self.producer, d) for d in producerargs ) )
        
    def debugwrap(self, func, dat):
        print >> sys.stderr, ".",
        d = func(dat)
        print >> sys.stderr, "+",
        return d
