
augroup coc_dynamic_autocmd
  autocmd!
  autocmd BufReadCmd,FileReadCmd,SourceCmd output:/* call coc#rpc#request('CocAutocmd', ['BufReadCmd','output', expand('<amatch>')])
  autocmd BufReadCmd,FileReadCmd,SourceCmd coc:/* call coc#rpc#request('CocAutocmd', ['BufReadCmd','coc', expand('<amatch>')])
  autocmd OptionSet runtimepath call coc#rpc#notify('OptionSet',[expand('<amatch>'), v:option_old, v:option_new])
  autocmd OptionSet completeopt call coc#rpc#notify('OptionSet',[expand('<amatch>'), v:option_old, v:option_new])
augroup end