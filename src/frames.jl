module Frames

using ..SimCore
using Topology


export frame_mechanism, write_popped_data_to_processor_memory, step_processor_one_cycle, read_data_from_processor_memory


frame_mechanism(framemover::Nothing, status, i) = return


"""
    frame_mechanism(framemover::FrameMover, status, i)

Simulates one cycle of the data plane for a given node `i`.

This function orchestrates the movement and processing of data frames for a single
localtick. It performs three main actions in sequence:
1.  **Receive**: It inspects the FIFOs for all incoming edges of node `i`.
    - It discards any frames that are older than the UGN time.
    - If a frame's arrival time (`send_tick + ugn`) matches the current
      `localtick`, it is popped from its FIFO and its data is written to the
      node's processor memory.
2.  **Process**: It signals the node's `framedataprocessor` (e.g., an Elara
    core model) to execute a single computation cycle.
3.  **Send**: It queries the `framedataprocessor` for data to be sent on each
    outgoing edge. A new `Frame` is created with this data, timestamped with
    the current `localtick`, and pushed into the corresponding outgoing FIFO.

# Arguments
- `framemover::FrameMover`: A struct containing the network's FIFOs.
- `status`: The overall simulation status, providing access to the model
   graph, UGN values, processor, and the node's local time (`theta`).
- `i`: The index of the node to be processed.
"""
function frame_mechanism(framemover::FrameMover, status, i)
    fdp = status.model.framedataprocessor
    graph = status.model.graph
    fifos = framemover.fifos
    ugns = status.model.ugns
    localtick = UInt64(floor(status.theta[i].y[end]))

    # get data out of Fifo and give it to fdp (ie Elara)
    for e in graph.incoming_edges[i]
        check_and_remove_old_frames(fifos[e], ugns[e], localtick)
        if first_frame_arrived(fifos[e], ugns[e], localtick)
            data::UInt64 = popfirst!(fifos[e]).data
            port = inportid_from_edgeid(graph, i, e)
            # send data to Elara processor
            write_popped_data_to_processor_memory(fdp, i, port, localtick, data)
        end
    end

    # do compute step
    step_processor_one_cycle(fdp, i, localtick)

    # get data from compute and put it onto the link
    for e in graph.outgoing_edges[i]
        port = outportid_from_edgeid(graph, i, e)
        # call Elara to get data from core/memory
        data::UInt64 = read_data_from_processor_memory(fdp, i, port, localtick)
        push!(framemover.fifos[e], Frame(localtick, data))
        # append to end of fifo
    end
end

function remove_old_frames(fifo, ugn, localtick)
    j = findlast(a ->  a[1] + ugn <= localtick, fifo)
    for k=1:j-1
        popfirst!(fifo)
    end
end

function check_and_remove_old_frames(fifo, ugn, localtick)
    if length(fifo) == 0
        return
    end
    if fifo[1].send_tick + ugn < localtick
        remove_old_frames(fifo, ugn, localtick)
    end
end

function first_frame_arrived(fifo, ugn, localtick)
    if isempty(fifo)
        return false
    end
    if fifo[1].send_tick + ugn != localtick
        return false
    end
    return true
end

write_popped_data_to_processor_memory(fdp::Nothing, i, port, localtick, data) = return
step_processor_one_cycle(fdp::Nothing, i, localtick) = return
read_data_from_processor_memory(fdp::Nothing, i, port, localtick)::UInt64 = 0x0



end
