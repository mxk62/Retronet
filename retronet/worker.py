import json
import zmq
from retronet import Chemical


def worker(num, transforms):
    """
    """

    context = zmq.Context()

    # Set up channel to request for a new task.
    receiver = context.socket(zmq.PULL)
    receiver.connect('tcp://localhost:5555')

    # Set up a channel to send the results to.
    sender = context.socket(zmq.PUSH)
    sender.connect('tcp://localhost:5556')

    # Set up channel for control input.
    controller = context.socket(zmq.SUB)
    controller.connect('tcp://localhost:5557')
    controller.setsockopt_string(zmq.SUBSCRIBE, ''.decode('ascii'))

    # Set up poller to multiplex the task and control receiver.
    poller = zmq.Poller()
    poller.register(receiver, zmq.POLLIN)
    poller.register(controller, zmq.POLLIN)

    while True:
        socks = dict(poller.poll())

        if socks.get(receiver) == zmq.POLLIN:
            # Request for a task.
            #request = 'FREE'
            #receiver.send_string(request)
            #print 'Sending request \'%s\'' % request

            # Process the task, i.e. make a retrostep.
            message = json.loads(receiver.recv_json())
            smiles, trans_id = message.values()

            chem = Chemical(smiles)
            print 'Worker-%d: applying transform %d to %s' % (
                num, trans_id, chem.smiles)
            reactant_sets = chem.make_retrostep(transforms[trans_id])

            # Send back the results.
            print 'Worker-%d: sending back %d results.' % (
                num, len(reactant_sets))
            message = {'results': reactant_sets if reactant_sets else 'NONE'}
            sender.send_json(json.dumps(message))

        if socks.get(controller) == zmq.POLLIN:
            message = controller.recv_string()
            if message == 'DONE':
                break