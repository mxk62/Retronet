import zmq


def worker(num):
    """
    """

    print 'Spawning %d worker' % num
    context = zmq.Context()

    # Set up channel to request for a new task.
    receiver = context.socket(zmq.PULL)
    receiver.connect('tcp://localhost:5555')
    print 'Receiver set up.'

    # Set up a channel to send the results to.
    sender = context.socket(zmq.PUSH)
    sender.connect('tcp://localhost:5556')
    print 'Sender set up.'

    # Set up channel for control input.
    controller = context.socket(zmq.SUB)
    controller.connect('tcp://localhost:5557')
    controller.setsockopt_string(zmq.SUBSCRIBE, ''.decode('ascii'))

    # Set up poller to multiplex the task and control receiver.
    poller = zmq.Poller()
    poller.register(receiver, zmq.POLLIN)
    poller.register(controller, zmq.POLLIN)
    print 'Poller set up.'

    while True:
        socks = dict(poller.poll())

        if socks.get(receiver) == zmq.POLLIN:
            # Request for a task.
            #request = 'FREE'
            #receiver.send_string(request)
            #print 'Sending request \'%s\'' % request

            # Process the task, i.e. make a retrostep.
            message = receiver.recv_pyobj()
            chem, transform = message
            print 'Received task: apply %s to %s' % (
                transform.smarts, chem.smiles)
            reactant_sets = chem.make_retrostep(transform)

            # Send back the results.
            print 'Worker %d: sending back %d results.' % (
                num, len(reactant_sets))
            sender.send_pyobj(reactant_sets)

        if socks.get(controller) == zmq.POLLIN:
            message = controller.recv_string()
            if message == 'DONE':
                break