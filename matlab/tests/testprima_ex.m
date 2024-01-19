function testprima_ex()
%TESTPRIMA tests prima extensively on a few VERY simple problems, combining testprima and pdv.

ver;

root_dir = fileparts(fileparts(pwd()));
cd(root_dir);
options = struct();
options.debug=true;
setup(options);
prima('info')
testprima(false, 1.0e-10, 100);
setup
setup path
prima('info')
testprima(false, 1.0e-10, 100);
setup cobyla
setup uobyqa
setup newuoa
setup bobyqa
setup lincoa
setup prima
setup path
setup clean
setup path
setup uninstall
setup path
setup uninstall
cd(fullfile(root_dir, 'matlab', 'tests'));
pdv
cd(root_dir);
setup
prima('info')
cd(fullfile(root_dir, 'matlab', 'examples'));
rosenbrock_example

end
