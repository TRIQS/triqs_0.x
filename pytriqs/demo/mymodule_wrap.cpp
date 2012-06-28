#include <boost/python.hpp>
#include <iostream>

using namespace boost::python;

class myClass
{
  public:
    myClass(object ob){ U = extract<int>(ob.attr("U")); }
    void solve() { std::cout << "Je suis dans le C++\n" << "U = " << U << std::endl; }
    int U;
};

BOOST_PYTHON_MODULE(mymodule)
{
    class_<myClass>("myClass", init<object>())
      .def("solvecpp",&myClass::solve)
    ;
};
