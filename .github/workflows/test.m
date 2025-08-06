options.blacklist = {};
  switch 'quadruple'
  case 'single'
      switch 'newuoa'
      case 'newuoa'
      case 'bobyqa'
      case 'lincoa'
          options.blacklist = [options.blacklist, {'EQC', 'GMNCASE3', 'PENTAGON'}];
      end
  case 'quadruple'
      switch 'newuoa'
      case 'newuoa'
          options.blacklist = [options.blacklist, {'ERRINROS', 'TOINTPSP', 'ARGLINC', 'LUKSAN13LS'}];
      case 'bobyqa'
          options.blacklist = [options.blacklist, {'ERRINRSM', 'TOINTGOR', 'DECONVU'}];
      case 'lincoa'
          options.blacklist = [options.blacklist, {'ERRINRSM', 'METHANB8LS', 'GMNCASE2'}];
      end
  end
