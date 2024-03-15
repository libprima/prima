function isms = ismac_silicon()
% ISMAC_SILICON checks whether the current machine is an Apple Silicon Mac.
isms = false;
if ismac
    [~, result] = system('uname -v');
    isms = contains(result, 'arm64', 'IgnoreCase', true);
end
