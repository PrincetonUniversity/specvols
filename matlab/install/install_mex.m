% INSTALL_MEX Compile MEX file part of the package
%
% Usage
%    install_mex();
%
% Description
%    Compiles all MEX files that are part of the package. This function should
%    only be called once when the package is first installed.

function install_mex()
    mex_dir = fullfile(pkg_root(), 'matlab', 'utils', 'mex');

    mex_files = dir(fullfile(mex_dir, '*.c'));

    current_dir = pwd();

    cd(mex_dir);

    if ~isoctave()
        mex_opts = {'-lgomp'};
    else
        mex_opts = {'-lgomp'};
    end

    for k = 1:numel(mex_files)
        filename =  mex_files(k).name;

        fprintf('Compiling ''%s''...\n', filename);
        try
            mex(mex_opts{:}, fullfile(mex_dir, filename));
        catch
            err_msg = sprintf('Unable to compile ''%s''.', filename);
            warning(err_msg);
        end
    end

    cd(current_dir);
end
