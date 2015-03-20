from multiprocessing import Pool


class ProducerConsumer(object):
    def __init__(self, args, threads, producer, consumer):
        self.threads= threads
        self.producer = producer
        self.consumer = consumer

    def run(self, passthrough, producerargs):
        pool = Pool(processes = threads)
        ret = self.pool.imap_unordered(self.producer, producerargs)
        self.consumer(passthrough, ret)
        pool.close()

        
