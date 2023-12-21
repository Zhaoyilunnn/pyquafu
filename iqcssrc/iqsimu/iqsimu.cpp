#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "simulator.hpp"
#include "mps_simulator.hpp"
#include "optimize/reduce.hpp"
#ifdef _USE_GPU
#include <cuda_simulator.cuh>
#endif

#ifdef _USE_CUQUANTUM
#include <custate_simu.cuh>
#endif

namespace py = pybind11;
using namespace pybind11::literals;

template <typename T>
py::array_t<T> to_numpy(const std::tuple<T*, size_t> &src) {
    auto src_ptr = std::get<0>(src);
    auto src_size = std::get<1>(src);

    auto capsule = py::capsule(src_ptr, [](void* p) {
        delete [] reinterpret_cast<T*>(p);
    });
    return py::array_t<T>(
        src_size,
        src_ptr,
        capsule
    );
}

template <typename T, typename allocator>
py::array_t<T> to_numpy(std::vector<T, allocator> &&src) {
  vector<T, allocator>* src_ptr = new std::vector<T, allocator>(std::move(src));
  auto capsule = py::capsule(src_ptr, [](void* p) { delete reinterpret_cast<std::vector<T, allocator>*>(p); });
  return py::array_t<T>(
    src_ptr->size(),  // shape of array
    src_ptr->data(),  // c-style contiguous strides for vector
    capsule           // numpy array references this parent
  );
}

py::object applyop_statevec(py::object const& pyop, py::array_t<complex<double>> &np_inputstate){
    py::buffer_info buf = np_inputstate.request();
    auto* data_ptr = reinterpret_cast<std::complex<double>*>(buf.ptr);
    size_t data_size = buf.size;

    auto op =  from_pyops(pyop);
    if (data_size == 0){
        return np_inputstate;
    }
    else{
        StateVector<double> state(data_ptr, buf.size);
        apply_op(state, op);
        state.move_data_to_python();
        return np_inputstate;
    }
}

py::object simulate_statevec(py::object const&pycircuit, int shots, py::array_t<complex<double>> &np_inputstate, bool merge_gates, int max_merge_size){
    auto get_full_mat = merge_gates;
    auto circuit = Circuit(pycircuit, get_full_mat, true);
    py::buffer_info buf = np_inputstate.request();
    auto* data_ptr = reinterpret_cast<std::complex<double>*>(buf.ptr);
    size_t data_size = buf.size;
    
    if (merge_gates){
        circuit = reduce_merge(circuit, max_merge_size);
    }

    if (data_size == 0){
        StateVector<double> state;
        auto counts = simulate(circuit, state, shots);
        py::dict pyres("counts"_a=counts, "statevector"_a=to_numpy(state.move_data_to_python()));
        // return to_numpy(state.move_data_to_python());
        return pyres;
    }
    else{
      StateVector<double> state(data_ptr, buf.size);
      auto counts = simulate(circuit, state, shots);
      state.move_data_to_python();
      py::dict pyres("counts"_a=counts, "statevector"_a= np_inputstate);
      return pyres;
    }
}


py::object simulate_mps(py::object const&pycircuit, py::list const paulis=py::list(0), double cutoff=1E-16,  int shots=0, bool save_statevec=false){
    auto circuit = Circuit(pycircuit, true, false);
    std::vector<std::map<int, char>> paulis_maps;
    for (auto pauli_h : paulis){
        py::object pypauli = py::reinterpret_borrow<py::object>(pauli_h);
        std::vector<pos_t> posv = pypauli.attr("pos").cast<std::vector<pos_t>>();
        string paulistr = pypauli.attr("paulistr").cast<string>();
        std::map<int, char> paulis_map;
        for (auto i = 0; i < posv.size();++i){
            paulis_map[posv[i]] = paulistr[i];
        }
        paulis_maps.push_back(paulis_map);
    }

    auto res = IQMPS::simulate_mps(circuit, paulis_maps, cutoff, shots, save_statevec);
    py::dict pyres("counts"_a=res.counts, "statevector"_a=to_numpy(std::move(res.ITVec)), "max_truncerr"_a=res.max_truncerr, "pauli_expects"_a=res.paulis_expects);
    return pyres;
}

py::object expect_statevec(py::array_t<complex<double>> const&np_inputstate, py::list const paulis)
{   
    py::buffer_info buf = np_inputstate.request();
    auto* data_ptr = reinterpret_cast<std::complex<double>*>(buf.ptr);
    size_t data_size = buf.size;
    StateVector<double> state(data_ptr, buf.size);
    py::list pyres;
    for (auto pauli_h : paulis){
         py::object pypauli = py::reinterpret_borrow<py::object>(pauli_h);
         std::vector<pos_t> posv = pypauli.attr("pos").cast<std::vector<pos_t>>();
         string paulistr = pypauli.attr("paulistr").cast<string>();
        double expec = state.expect_pauli(paulistr, posv);
        pyres.attr("append")(expec);
    }
    state.move_data_to_python();
    return pyres;
}

#ifdef _USE_GPU
py::object simulate_statevec_gpu(py::object const&pycircuit, py::array_t<complex<double>> &np_inputstate){
    auto circuit = Circuit(pycircuit);
    py::buffer_info buf = np_inputstate.request();
    auto* data_ptr = reinterpret_cast<std::complex<double>*>(buf.ptr);
    size_t data_size = buf.size;


    if (data_size == 0){
        StateVector<double> state;
        simulate_gpu(circuit, state);
        return to_numpy(state.move_data_to_python());
    }
    else{
      StateVector<double> state(data_ptr, buf.size);
      simulate_gpu(circuit, state);
      state.move_data_to_python();
      return np_inputstate;
    }
}
#endif

#ifdef _USE_CUQUANTUM
py::object simulate_statevec_custate(py::object const&pycircuit, py::array_t<complex<double>> &np_inputstate){
    auto circuit = Circuit(pycircuit);
    py::buffer_info buf = np_inputstate.request();
    auto* data_ptr = reinterpret_cast<std::complex<double>*>(buf.ptr);
    size_t data_size = buf.size;


    if (data_size == 0){
        StateVector<double> state;
        simulate_custate(circuit, state);
        return to_numpy(state.move_data_to_python());
    }
    else{
      StateVector<double> state(data_ptr, buf.size);
      simulate_custate(circuit, state);
      state.move_data_to_python();
      return np_inputstate;
    }
}
#endif


PYBIND11_MODULE(iqsimu, m) {
    m.doc() = "iqcs C++ simulator";
    
    m.def("simulate_statevec", &simulate_statevec, "Simulate with circuit", py::arg("circuit"), py::arg("shots") = 0, py::arg("inputstate")= py::array_t<complex<double>>(0), py::arg("merge_gates")=true, py::arg("max_merge_size")=5);

    m.def("simulate_mps", &simulate_mps, "Simulate circuit using MPS", py::arg("circuit"),  py::arg("paulis")=py::list(0), py::arg("cutoff")=1E-10,  py::arg("shots")=0, py::arg("save_statevec")=false);

    m.def("expect_statevec", &expect_statevec, "Calculate paulis expectation", py::arg("inputstate"), py::arg("paulis"));

    m.def("applyop_statevec", &applyop_statevec, "Apply single operator to state", py::arg("operation"), py::arg("inputstate"));

    #ifdef _USE_GPU
     m.def("simulate_statevec_gpu", &simulate_statevec_gpu, "Simulate with circuit", py::arg("circuit"), py::arg("inputstate")= py::array_t<complex<double>>(0));
    #endif

    #ifdef _USE_CUQUANTUM
    m.def("simulate_statevec_custate", &simulate_statevec_custate, "Simulate with circuit", py::arg("circuit"), py::arg("inputstate")= py::array_t<complex<double>>(0));
    #endif
}

