module OpenVerticalChannel

import EarthBox.ModelDataContainer: ModelData
import ..VerticalChannelBC

function set_open_vertical_channel!(model::ModelData)::Nothing
    VerticalChannelBC.set_open_vertical_channel!(model)
    return nothing
end

end # module OpenVerticalChannel 