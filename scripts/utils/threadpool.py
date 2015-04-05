from multiprocessing import Pool


class ProducerConsumer(object):
    def __init__(self, args, threads, producer, consumer):
        if threads > 1:
            self.threads= threads - 1
        else:
            self.threads = 1
        self.producer = producer
        self.consumer = consumer

    def run(self, passthrough, producerargs):
        pool = Pool(processes = self.threads)
        ret = pool.imap_unordered(self.producer, producerargs)
        self.consumer(passthrough, ret)
        pool.close()

        
