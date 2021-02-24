
#include <pybind11/stl.h>
namespace py = pybind11;
typedef std::map<int, double> dict;



dict
convert_dict_to_map_d(py::dict dictionary)
{
    std::map<int, double> result;
    for (std::pair<py::handle, py::handle> item : dictionary)
    {
        auto key = item.first.cast<int>();
        auto value = item.second.cast<double>();
        //cout << key << " : " << value;
        result[key] = value;
    }
    return result;
}
