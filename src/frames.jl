module Frames

using ..SimCore
using Topology


export frame_mechanism, received_data, process_step, data_to_send


frame_mechanism(framemover::Nothing, status, i) = return

# Function calls:
#  Controller -> vmtick -> receive_data -> checkfifo -> received_data
#                       -> send_data -> data_to_send
#
# Dataflow:
#  data_to_send -> send_data -> vmtick -> receive_data -> checkfifo -> received_data
#
#

#  send_data gets data from data_to_send and puts it in a frame on the fifo
#  data_to_send(fs, i, port, localtick) returns a UInt64
#
#  check_fifo gets a frame out of the fifo and gives the data field to received_data
#  received_data(fs, i, port, localtick, data::UInt64)  returns nothing
#
#  So Callisto needs a FIFO. The frame type can be (localtick, data)
#  It's optional in C1, and built in to C2.
#

# TODO: refactor this to bring the loops into this function
# and to separate the checkfifo function  from the received_data function
# so that checkfifo returns a status and we call received_data
# from here.
function frame_mechanism(framemover::FrameMover, status, i)
    fdp = status.model.framedataprocessor
    localtick = UInt64(floor(status.theta[i].y[end]))

    # get data out of Fifo and give it to fdp (ie Elara)
    receive_data(framemover, status, i, localtick)

    # do compute step
    process_step(fdp, i, localtick)

    # get data from compute and put it onto the link
    send_data(framemover, status, i, localtick)
end

function clearfifo(fifo, ugn, localtick)
    j = findlast(a ->  a[1] + ugn <= localtick, fifo)
    for k=1:j-1
        popfirst!(fifo)
    end
end

function receive_data(framemover::FrameMover, status, i, localtick)
    fdp = status.model.framedataprocessor
    graph = status.model.graph
    for e in graph.incoming_edges[i]
        port = inportid_from_edgeid(graph, i, e)
        fifo = framemover.fifos[e]
        if length(fifo) > 0
            ugn = status.model.ugns[e]
            #
            # checkfifo will pop the top frame from the fifo
            # and call received_data to handle it
            checkfifo(fdp, fifo, ugn, i, port, localtick)
        end
    end
end


function checkfifo(fdp, fifo, ugn, i, port, localtick)
    # remove frames that have arrived before now
    if fifo[1].send_tick + ugn < localtick
        clearfifo(fifo, ugn, localtick)
    end
    if isempty(fifo)
        return
    end
    # get the current frame and receive it's data
    if fifo[1].send_tick + ugn == localtick
        data = popfirst!(fifo).data
        received_data(fdp,i, port, localtick, data)
    end
end

function send_data(framemover::FrameMover, status, i, localtick)
    fdp = status.model.framedataprocessor
    graph = status.model.graph
    for e in graph.outgoing_edges[i]
        port = outportid_from_edgeid(graph, i, e)
        #
        # data_to_send is another Elara entrypoint.
        # It gets the data from memory so that Callisto
        # can put it into the fifo. Here data is a UInt64.
        #
        data = data_to_send(fdp, i, port, localtick)
        push!(framemover.fifos[e], Frame(localtick, data))  # append to end of fifo
    end
end

received_data(fdp::Nothing, i, port, localtick, data) = return
process_step(fdp::Nothing, i, localtick) = return
data_to_send(fdp::Nothing, i, port, localtick) = return







end

