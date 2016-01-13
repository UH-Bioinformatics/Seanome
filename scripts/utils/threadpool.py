from multiprocessing import Pool


class ProducerConsumer(object):
    def __init__(self, args, threads, producer, consumer):
        # we subtract 1 thread, to account for the consumer thread. aka the main thread
		self.debug = args.debug
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
			self.consumer(passthrough, ( self.producer(d) for d in producerargs ) )
				

        
