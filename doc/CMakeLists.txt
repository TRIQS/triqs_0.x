#option(Build_Documentation "Build the documentation ?" OFF)

option(Sphinx_Math_Generator_MathJax " " OFF)

if (Sphinx_Math_Generator_MathJax)
 SET(SPHINX_MATH_GENERATOR "mathjax")
else (Sphinx_Math_Generator_MathJax)
 SET(SPHINX_MATH_GENERATOR "pngmath")
endif (Sphinx_Math_Generator_MathJax)

add_subdirectory(user_manual)
add_subdirectory(developer_manual)

file(COPY ../foreignlibs/breathe ../foreignlibs/doxylink ../foreignlibs/mpl DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


