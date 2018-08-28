function GetNDVI,infile,ndvifile,error=error
  compile_opt idl2
  
  call_procedure,'envi_open_file',infile,r_fid=tfid
  if tfid eq -1 then begin
    error=infile+'文件打开失败!'
    return,0
  endif
  call_procedure,'envi_file_query',tfid,dims=dims,ns=ns,nl=nl,nb=nb
  mapinfo=call_function('envi_get_map_info',fid=tfid)
  
  expstr='(b4*1.0-b3*1.0)/(b4*1.0+b3*1.0)'
  tempfile=file_dirname(infile)+file_basename(ndvifile,'.tiff')+'_.dat'
  call_procedure,'envi_doit','math_doit', $
     fid=[tfid,tfid], $
     dims=dims, $
     pos=[2,3], $
     exp=expstr, $
     out_bname='NDVI', $
     out_name=tempfile, $
     r_fid=rfid  
  if rfid eq -1 then begin
    error='NDVI计算失败!'
    return,0
  endif
  call_procedure,'envi_file_mng',id=tfid,/remove
  call_procedure,'envi_file_mng',id=rfid,/remove
  
  call_procedure,'envi_open_file',tempfile,r_fid=tfid
  call_procedure,'envi_file_query',tfid,dims=dims,ns=ns,nl=nl,nb=nb,bnames=bnames
  call_procedure,'envi_output_to_external_format', $
    fid=tfid, $
    dims=dims, $
    pos=indgen(nb), $
    out_bname=bnames, $
    out_name=ndvifile, $
    /tiff
  call_procedure,'envi_file_mng', id=tfid, /remove, /delete
  if ~file_test(ndvifile) then begin
    error='NDVI结果转TIFF格式失败!'
    return,0
  endif  
  
  return,1  
end