function Preprocess,orignalfile,outfile,error=error
  compile_opt idl2
;  envi,/restore_base_save_files
;  envi_batch_init
;  envi_batch_status_window,/on  
  e=call_function('envi',/headless)
  
  xmlfile=file_dirname(orignalfile)+path_sep()+file_basename(orignalfile,'.tiff')+'.xml'
  rpbfile=file_dirname(orignalfile)+path_sep()+file_basename(orignalfile,'.tiff')+'.rpb'
  if ~file_test(xmlfile) then begin
    error=xmlfile+'文件不存在!'
    return,0
  endif
  if ~file_test(rpbfile) then begin
    error=rpbfile+'文件不存在!'
    return,0
  endif
  
  ;读取原始数据相关信息
  xmlfile=file_dirname(orignalfile)+path_sep()+file_basename(orignalfile,'.tiff')+'.xml'
  ret=ReadInfo(xmlfile,satellite,sensor,year,month,day,gmt,latitude,longitude,latrange,lonrange,error=error)
  if ret ne 1 then begin
    error='读取原始数据相关信息失败!---'+error
    print,error
    return,0
  endif
  print,'Satellite:',satellite
  print,'Sensor:',sensor
  print,'Year:',year,'Month:',month,'Day:',day
  print,'Longitude:',longitude
  print,'Latitude:',latitude
  print,'Lonrange:',lonrange
  print,'Latrange:',latrange
  
  currentdir=file_dirname(file_dirname(routine_filepath('Preprocess',/is_function)))+path_sep()
  resourcedir=currentdir+'Resource'+path_sep()
  if ~file_test(resourcedir,/directory) then begin
    error=resourcedir+'辅助文件路径不存在!'
    return,0
  endif
  defsysv,'!resourcedir',resourcedir
  tempdir=currentdir+'Temp'+path_sep()
  if ~file_test(tempdir) then file_mkdir,tempdir
  defsysv,'!tempdir',tempdir
    
  ;读取定标系数  
  paramfile=currentdir+'Resource'+path_sep()+'Calibration_Parametres.xml'  
  ret=ReadParam(paramfile,year,satellite,sensor,gain,offset,esum,error=error)
  if ret ne 1 then begin
    error='读取定标系数失败!---'+error
    print,error
    return,0
  endif
  print,'Gain:',strjoin(strtrim(string(gain),2),',')
  print,'Offset:',strjoin(strtrim(string(offset),2),',')
  print,'Esum:',esum
    
  ;辐射定标
  radfile=!tempdir+file_basename(outfile,'.tiff')+'_Rad.dat'
  ret=Calibrate(year,satellite,sensor,orignalfile,gain,offset,radfile,error=error)  
  if ret ne 1 then begin
    error='辐射定标失败!---'+error
    print,error
    return,0
  endif
    
  ;Flaash大气校正  
  reffile=!tempdir+file_basename(outfile,'.tiff')+'_Reflectance.tiff'
  ret=Flaash(radfile,satellite,sensor,year,month,day,gmt,latitude,longitude,latrange,lonrange,reffile,error=error)
  if ret ne 1 then begin
    error='Flaash大气校正失败!---'+error
    print,error
    return,0
  endif 

  ret=Orthorpc(reffile,rpbfile,satellite,sensor,outfile,error=error)
  if ret ne 1 then begin
    error='正射校正失败!---'+error
    print,error
    return,0
  endif
  
  return,1
end


function ReadInfo,xmlfile,satellite,sensor,year,month,day,gmt,latitude,longitude,latrange,lonrange,error=error
  compile_opt idl2
  
  catch,error_status
  if error_status ne 0 then begin
    error=!ERROR_STATE.MSG
    catch,/cancel
    return,0
  endif

  if ~file_test(xmlfile) then begin
    error=xmlfile+'文件不存在!'
    return,0
  endif

  oDoc=obj_new('IDLffXMLDOMDocument', filename=xmlfile)

  oXmlEle=oDoc->GetDocumentElement()                          ;根节点
  
  oSatellite=oXmlEle->GetElementsByTagName('SatelliteID')
  Satellite_list=oSatellite->item(0)
  satellite=(Satellite_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oSatellite, Satellite_list]
  
  oSensor=oXmlEle->GetElementsByTagName('SensorID')
  Sensor_list=oSensor->item(0)
  sensor=(Sensor_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oSensor, Sensor_list]

  oInt=oXmlEle->GetElementsByTagName('ReceiveTime')
  Int_list=oInt->item(0)
  ReceiveTime=(Int_list->GetFirstChild())->GetNodeValue()
  date_=(strsplit(ReceiveTime,' ',/extract))[0]
  temp=fix(strsplit(date_,'-',/extract))
  year=temp[0] & month=temp[1] & day=temp[2]
  time=(strsplit(ReceiveTime,' ',/extract))[1]
  temp=strsplit(time,':',/extract)
  hour=temp[0] & minute=temp[1] & second=temp[2]
  gmt=float(hour)+float(minute)/60+float(second)/3600
  obj_destroy,[oInt, Int_list]

  oInt=oXmlEle->GetElementsByTagName('TopLeftLatitude')
  Int_list=oInt->item(0)
  TopLeftLatitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  oInt=oXmlEle->GetElementsByTagName('TopRightLatitude')
  Int_list=oInt->item(0)
  TopRightLatitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  oInt=oXmlEle->GetElementsByTagName('BottomRightLatitude')
  Int_list=oInt->item(0)
  BottomRightLatitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  oInt=oXmlEle->GetElementsByTagName('BottomLeftLatitude')
  Int_list=oInt->item(0)
  BottomLeftLatitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  latitude=(float(TopLeftLatitude)+float(TopRightLatitude)+float(BottomRightLatitude)+float(BottomLeftLatitude))/4
  latrange=[min([BottomRightLatitude,BottomLeftLatitude]),max([TopLeftLatitude,TopRightLatitude])]

  oInt=oXmlEle->GetElementsByTagName('TopLeftLongitude')
  Int_list=oInt->item(0)
  TopLeftLongitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  oInt=oXmlEle->GetElementsByTagName('TopRightLongitude')
  Int_list=oInt->item(0)
  TopRightLongitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  oInt=oXmlEle->GetElementsByTagName('BottomRightLongitude')
  Int_list=oInt->item(0)
  BottomRightLongitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  oInt=oXmlEle->GetElementsByTagName('BottomLeftLongitude')
  Int_list=oInt->item(0)
  BottomLeftLongitude=(Int_list->GetFirstChild())->GetNodeValue()
  obj_destroy,[oInt, Int_list]
  longitude=(float(TopLeftLongitude)+float(TopRightLongitude)+float(BottomRightLongitude)+float(BottomLeftLongitude))/4
  lonrange=[min([TopLeftLongitude,BottomLeftLongitude]),max([TopRightLongitude,BottomRightLongitude])]

  obj_destroy, oDoc  
  return,1
end


function ReadParam,parametresfile,year,satellite,sensor,gain,offset,Esum,error=error
  compile_opt idl2

  catch,error_status
  if error_status ne 0 then begin
    error=!ERROR_STATE.MSG
    catch,/cancel
    return,0
  endif
    
  if ~file_test(parametresfile) then begin
    error=parametresfile+'参数文件不存在!'
    return,0
  endif

  oDoc=obj_new('IDLffXMLDOMDocument', filename=parametresfile)
  oxmlele=oDoc->GetDocumentElement()                                ;根节点
  oyear=oXmlEle->GetElementsByTagName('Year')
  num0=oyear->getlength()
  for i=0, num0-1 do begin                                          ;年份
    
    year_=oyear->item(i)
    if year eq year_->GetAttribute('value') then begin
      osatellite=year_->GetElementsByTagName('Satellite')
      num1=osatellite->getlength()
      for j=0, num1-1 do begin                                      ;卫星
        
        satellite_=osatellite->item(j)
        if satellite eq satellite_->GetAttribute('name') then begin
          osensor=satellite_->GetElementsByTagName('Sensor')
          num2=osensor->getlength()
          for k=0, num2-1 do begin                                  ;传感器
            
            sensor_=osensor->item(k)
            if sensor eq sensor_->GetAttribute('name') then begin
              oband=sensor_->GetElementsByTagName('Element')
              num3=oband->getlength()
              gain=fltarr(num3)
              offset=fltarr(num3)
              Esum=fltarr(num3)
              for t=0, num3-1 do begin                               ;波段
                band_=oband->item(t)
                value=band_->GetAttribute('value')
                gain[t]=float((strsplit(value,',',/extract))[0])
                offset[t]=float((strsplit(value,',',/extract))[1])
                Esum[t]=float((strsplit(value,',',/extract))[2])
              endfor
              obj_destroy, oDoc
              break
            endif
            
          endfor
          break
        endif
        
      endfor
      break
    endif
        
  endfor
  
  if (n_elements(gain) eq 0) or (n_elements(offset) eq 0) or (n_elements(esum) eq 0) then begin
    error='没有符合要求的定标系数,读取失败!'
    return,0
  endif

  if obj_valid(oDoc) then obj_destroy, oDoc  
  return,1
end


function Calibrate,year,satellite,sensor,inputfile,gain,offset,radfile,error=error
  compile_opt idl2  

  catch,error_status
  if error_status ne 0 then begin
    error=!ERROR_STATE.MSG
    catch,/cancel
    return,0
  endif  

  if ~file_test(inputfile) then begin
    error='原始数据不存在!
    return,0
  endif
  if n_elements(gain) eq 0 then begin
    error='增益值未指定!'
    return,0
  endif
  if n_elements(offset) eq 0 then begin
    error='偏差值未指定!'
    return,0
  endif

  ;辐射定标
  call_procedure, 'envi_open_file', inputfile, r_fid=tfid
  call_procedure, 'envi_file_query', tfid, dims=dims, ns=ns, nl=nl, nb=nb
  mapinfo0=call_function('envi_get_map_info',fid=tfid)
  oproj=call_function('envi_get_projection',fid=tfid)
  fidarr=!null & raddata=!null

  tempfiles=!null
  for i=0, nb-1 do begin
    ;不同年份,不同卫星,定标公式不同
    case satellite of
      'GF1':begin
        if year eq '2013' then begin
          if strmid(sensor,0,3) eq 'PMS' then begin
            expstr='b1*('+strtrim(string(gain[i]),2)+')+('+strtrim(string(offset[i]),2)+')'
          endif else if strmid(sensor,0,3) eq 'WFV' then begin
            expstr='(b1-('+strtrim(string(offset[i]),2)+'))/('+strtrim(string(gain[i]),2)+')'
          endif
        endif else if year eq '2014' then begin
          expstr='b1*('+strtrim(string(gain[i]),2)+')'
        endif else begin
          expstr='b1*('+strtrim(string(gain[i]),2)+')+('+strtrim(string(offset[i]),2)+')'
        endelse
      end
      'GF2':expstr='b1*('+strtrim(string(gain[i]),2)+')+('+strtrim(string(offset[i]),2)+')'
      else:
    endcase

    tempfile=!tempdir+file_basename(inputfile,'.tiff')+'_'+strtrim(string(i),2)+'.dat'
    call_procedure, 'envi_doit', 'math_doit', $
      exp=expstr, $
      fid=[tfid], $
      dims=dims, $
      pos=[i], $
      out_name=tempfile, $
      r_fid=rfid
    if rfid eq -1 then begin
      error=strtrim(string(i+1),2)+'波段辐射定标失败!'
      return,0
    endif
    tempfiles=[tempfiles,tempfile]
    call_procedure, 'envi_file_mng', id=rfid, /remove
  endfor

  radfile_=!tempdir+file_basename(radfile,'.dat')+'_.dat'
  openw,lun,radfile_,/get_lun

  for i=0, nb-1 do begin
    call_procedure, 'envi_open_file', tempfiles[i], r_fid=tempfid
    raddata=call_function('envi_get_data',fid=tempfid,dims=dims,pos=[0])
    writeu,lun,temporary(raddata)
    fidarr=[fidarr,tempfid]
  endfor
  free_lun,lun
  
  call_procedure, 'envi_file_query', tempfid, data_type=datatype
  call_procedure, 'envi_setup_head', $
    fname=radfile_, $
    ns=ns,nl=nl,nb=nb, $
    offset=0, $
    interleave=0, $
    data_type=datatype, $
;    map_info=mapinfo, $
    r_fid=rfid, $
    /write
  if rfid eq -1 then begin
    error='波段组合失败!'
    return,0
  endif
  call_procedure,'envi_file_mng',id=rfid,/remove
  call_procedure,'envi_open_file',radfile_,r_fid=rfid
  call_procedure,'envi_file_query',rfid,dims=dims,ns=ns,nl=nl,nb=nb
  call_procedure,'envi_doit','convert_doit', $
    fid=rfid, $
    dims=dims, $
    pos=indgen(nb), $
    o_interleave=1, $
    out_name=radfile, $
    r_fid=rfid_
  if rfid_ eq -1 then begin
    error='BSQ转BIL失败!'
    return,0
  endif

  for i=0, nb-1 do begin
    call_procedure, 'envi_file_mng', id=fidarr[i], /remove, /delete
  endfor
  call_procedure, 'envi_file_mng', id=tfid, /remove;, /delete
  call_procedure, 'envi_file_mng', id=rfid, /remove, /delete
  call_procedure, 'envi_file_mng', id=rfid_, /remove

  return,1  
end


function Ground_Elevation,demfile,latrange,lonrange,meandem,error=error
  compile_opt idl2
  
  call_procedure, 'envi_open_file', demfile, r_fid=tfid
  if tfid eq '-1' then begin
    error='获取目标区域海拔时DEM文件打开失败!'
    return,0
  endif
  call_procedure, 'envi_file_query', tfid, dims=dims, ns=ns, nl=nl, nb=nb
  iproj=call_function('envi_proj_create', /geographic)
  oproj=call_function('envi_get_projection', fid=tfid)
  call_procedure, 'envi_convert_projection_coordinates', $
    lonrange, latrange, iproj, $
    xmap, ymap, oproj
  call_procedure, 'envi_convert_file_coordinates', $
    tfid, xf, yf, xmap, ymap

  sub_dims=[-1, round(xf[0]), round(xf[1]), round(yf[1]), round(yf[0])]
  print, sub_dims
  tempfile=file_dirname(demfile)+path_sep()+'DEM_Subset.dat'
  call_procedure, 'envi_doit', 'resize_doit', $
    fid=tfid, $
    dims=sub_dims, $
    pos=indgen(nb), $
    interp=0, $
    out_name=tempfile, $
    rfact=[1.,1.], $
    r_fid=rfid
  if rfid eq -1 then begin
    error='获取目标区域平均海拔失败!'
    return,0
  endif
  
  call_procedure,'envi_file_mng',id=rfid,/remove
  call_procedure,'envi_open_file',tempfile,r_fid=rfid
  call_procedure, 'envi_file_query', rfid, dims=dims, ns=ns, nl=nl, nb=nb
  demdata=call_function('envi_get_data', fid=rfid, dims=dims, pos=[0])
  index=where(demdata ne 0, count)
  if count gt 0 then begin
    meandem=mean(demdata[index])/1000.
  endif else begin
    meandem=0.
  endelse
  demdata=!Null

  call_procedure, 'envi_file_mng', id=tfid, /remove
  call_procedure, 'envi_file_mng', id=rfid, /remove, /delete

  return,1  
end


function Atmosphere_Model,month,latitude,atmos_model,error=error
  compile_opt idl2
  
  ;SAW-0, MLW-1, US.Standard-2, SAS-3, MLS-4, T-5
  atmosphere_models= $
    ; 1&2    3&4    5&6    7&8    9&10   11&12
    [['SAW', 'SAW', 'SAW', 'MLW', 'MLW', 'SAW'], $   ;75~85
    ['SAW', 'SAW', 'MLW', 'MLW', 'MLW', 'SAW'], $   ;65~75
    ['MLW', 'MLW', 'MLW', 'SAS', 'SAS', 'MLW'], $   ;55~65
    ['MLW', 'MLW', 'SAS', 'SAS', 'SAS', 'SAS'], $   ;45~55
    ['SAS', 'SAS', 'SAS', 'MLS', 'MLS', 'SAS'], $   ;35~45
    ['MLS', 'MLS', 'MLS', 'T',   'T',   'MLS'], $   ;25~35
    ['T',   'T',   'T',   'T',   'T',   'T'  ], $   ;15~25
    ['T',   'T',   'T',   'T',   'T',   'T'  ], $   ;05~15
    ['T',   'T',   'T',   'T',   'T',   'T'  ], $   ;-05~05
    ['T',   'T',   'T',   'T',   'T',   'T'  ], $   ;-15~-05
    ['T',   'T',   'T',   'MLS', 'MLS', 'T'  ], $   ;-25~-15
    ['MLS', 'MLS', 'MLS', 'MLS', 'MLS', 'MLS'], $   ;-35~-25
    ['SAS', 'SAS', 'SAS', 'SAS', 'SAS', 'SAS'], $   ;-45~-35
    ['SAS', 'SAS', 'SAS', 'MLW', 'MLW', 'SAS'], $   ;-55~-45
    ['MLW', 'MLW', 'MLW', 'MLW', 'MLW', 'MLW'], $   ;-65~-55
    ['MLW', 'MLW', 'MLW', 'MLW', 'MLW', 'MLW'], $   ;-75~-65
    ['MLW', 'MLW', 'MLW', 'MLW', 'MLW', 'MLW'] ]

  if fix(month) mod 2 eq 1 then sam=(fix(month)+1)/2-1 else sam=fix(month)/2-1
  ;8 7 6 5 4 3 2 1 0 -1 -2 -3 -4 -5 -6 -7 -8
  ;0 1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16
  if latitude-fix(latitude/10)*10 ge 5 then lin=8-(fix(latitude)/10+1) else lin=8-(fix(latitude)/10)
  case atmosphere_models[sam, lin] of
    'SAW': atmos_model=0
    'MLW': atmos_model=1
    'US':  atmos_model=2
    'SAS': atmos_model=3
    'MLS': atmos_model=4
    'T':   atmos_model=5
  endcase
  
  return,1    
end


function Readfwhm,fwhmfile,satellite,sensor,wavelength,fwhm,error=error
  compile_opt idl2
  
  if ~file_test(fwhmfile) then begin
    error=fwhmfile+'文件不存在!'
    return,0
  endif
  
  oDoc=obj_new('IDLffXMLDOMDocument', filename=fwhmfile)
  oxmlele=oDoc->GetDocumentElement()                            ;根节点
  osatellite=oxmlele->GetElementsByTagName('Satellite')
  num1=osatellite->getlength()
  for j=0, num1-1 do begin                                      ;卫星
    satellite_=osatellite->item(j)
    if satellite eq satellite_->GetAttribute('name') then begin
      osensor=satellite_->GetElementsByTagName('Sensor')
      num2=osensor->getlength()
      for k=0, num2-1 do begin                                  ;传感器
        sensor_=osensor->item(k)
        if sensor eq sensor_->GetAttribute('name') then begin
          oband=sensor_->GetElementsByTagName('Element')
          num3=oband->getlength()
          wavelength=fltarr(num3)
          fwhm=fltarr(num3)
          for t=0, num3-1 do begin                               ;波段
            band_=oband->item(t)
            value=band_->GetAttribute('value')
            wavelength[t]=float((strsplit(value,',',/extract))[0])
            fwhm[t]=float((strsplit(value,',',/extract))[1])
          endfor
          obj_destroy, oDoc
          break
        endif
      endfor
      break
    endif
  endfor

  if obj_valid(oDoc) then obj_destroy, oDoc

  return,1   
end


function Flaash,radfile,satellite,sensor,year,month,day,gmt,latitude,longitude,latrange,lonrange,reffile,error=error
  compile_opt idl2

  catch,error_status
  if error_status ne 0 then begin
    error=!ERROR_STATE.MSG
    catch,/cancel
    return,0
  endif  
  
  call_procedure,'envi_open_file',radfile,r_fid=tfid
  if tfid eq -1 then begin
    error=radfile+'打开失败!'
    return,0
  endif
  call_procedure,'envi_file_query',tfid,dims=dims,ns=nspatial,nl=nlines,data_type=data_type

  reffile_=!tempdir+file_basename(reffile,'.tiff')+'_.dat'

  ;输出路径
  modtran_directory=file_dirname(reffile)+path_sep()

  ;光谱响应函数
  filter_func_filename=!resourcedir+'Filter'+path_sep()+satellite+'_'+sensor+'_'+'SpectralResponsivity.sli'
  if ~file_test(filter_func_filename) then begin
    error=filter_func_filename+'光谱响应函数不存在!'
    return,0
  endif

  ;读取XML文件获取影像时间等信息
;  ret=self->ObtainInfo(year,month,day,gmt,latitude,longitude,latrange,lonrange)

  ;卫星高度和影像分辨率
  case satellite of
    'GF1':begin
      if strmid(sensor,0,3) eq 'PMS' then pixel_size=8. else pixel_size=16.
      sensor_altitude=645.0
      sensor_name='GF-1'
      filter_func_file_index=0  ;????????????
    end
    'GF2':begin
      pixel_size=4.
      sensor_altitude=631.0
      sensor_name='GF-2'
      filter_func_file_index=0
    end
    else:begin
      error='satellite有误!'+satellite
      return,0
    endelse
  endcase

  ;研究区平均海拔高度
  demfile=!resourcedir+'GMTED2010.jp2'
  res=Ground_Elevation(demfile,latrange,lonrange,elevation,error=error)
  if res ne 1 then begin
    error='统计平均海拔高度失败!'+error
    return,0  
  endif
  ground_elevation=elevation

  ;大气模型：0-SAW;1-MLW;2-U.S. Standard;3-SAS;4-MLS;5-T
  res=Atmosphere_Model(month,latitude,atmosphere,error=error)
  if res ne 1 then begin
    error='大气模型选取失败!'+error
    return,0
  endif
  atmosphere_model=atmosphere

  ;气溶胶模型：0-No Aerosol;1-Rural;2-Maritime;3-Urban;4-Tropospheric
  aerosol_model=1

  ;波长相关信息
  wavelength_units='micron'
  fwhmfile=!resourcedir+'FWHM_Parametres.xml'
  if ~file_test(fwhmfile) then begin
    error=fwhmfile+'文件不存在!'
    return,0
  endif
  ret=Readfwhm(fwhmfile,satellite,sensor,wavelength,fwhm,error=error)
  if ret ne 1 then begin
    error='读取FWHM失败!'+error
    return,0
  endif

  ;辐射亮度单位转换系数
  input_scale=make_array(4,value=10.0,/double)

  flaash_obj=obj_new('flaash_batch', /anc_delete)

  ;设置大量的输入参数
  flaash_obj->SetProperty, $
    hyper = 0, $                     ;设置为1，表示高光谱；设置为0，表示多光谱
    radiance_file = radfile, $
    reflect_file = reffile_, $
    filter_func_filename = filter_func_filename, $
    filter_func_file_index = filter_func_file_index, $
    water_band_choice = 1.13, $
    red_channel = 3, $               ;0表示undefined，LC8红波段为第4波段
    green_channel = 2, $             ;0表示undefined，LC8绿波段为第3波段
    blue_channel = 0, $              ;0表示undefined，LC8蓝波段为第2波段
    water_retrieval = 0, $           ;Water Retrieval参数。0表示No，1表示Yes
    water_abs_channel = 0, $
    water_ref_channel = 0, $
    kt_upper_channel = 0, $          ;设置短波红外2（SWIR 2）
    kt_lower_channel = 3, $          ;设置红波段（Red）
    kt_cutoff = 0.08, $              ;Maximum Upper Channel Reflectance
    kt_ratio = 0.500, $              ;Reflectance Ratio
    cirrus_channel = 0, $            ;0表示undefined
    modtran_directory = modtran_directory, $
    visvalue = 40.0000, $            ;能见度，默认40km
    f_resolution = 5.0000, $
    day = day, $
    month = month, $
    year = year, $
    gmt = gmt, $
    latitude = latitude, $
    longitude = longitude, $
    sensor_altitude = sensor_altitude, $   ;传感器高度
    ground_elevation = ground_elevation, $ ;平均海拔，单位km
    view_zenith_angle = 180, $
    view_azimuth = 0, $
    atmosphere_model = atmosphere_model, $       ;大气模型：0-SAW;1-MLW;2-U.S. Standard;3-SAS;4-MLS;5-T
    aerosol_model = aerosol_model, $             ;气溶胶模型：0-No Aerosol;1-Rural;2-Maritime;3-Urban;4-Tropospheric
    multiscatter_model = 2, $
    disort_streams = 8, $
    co2mix = 390.0000, $
    water_column_multiplier = 1.0000, $
    nspatial = nspatial, $
    nlines = nlines, $
    data_type = data_type, $
    margin1 = 0, $
    margin2 = 0, $
    nskip = 0, $
    pixel_size = pixel_size, $
;    sensor_name=sensor_name, $
    sensor_name='UNKNOWN-MSI', $
    aerosol_scaleht = 1.5000, $
    use_adjacency = 1, $             ;中高分辨率设置为1，低分辨率（如Modis）设置为0
    output_scale = 10000.0000, $     ;输出结果缩放系数
    polishing_res = 0, $             ;对应 Width (number of bands) 参数，多光谱设置0即可。
    aerosol_retrieval = 0, $         ;0 表示 None；1 表示 2-Band (K-T)；2 表示 2-Band Over Water
    calc_wl_correction = 0, $        ;对应FLAASH面板中的 Wavelength Recalibration，多光谱一般为0
    reuse_modtran_calcs = 0, $
    use_square_slit_function = 0, $
    convolution_method = 'fft', $
    use_tiling = 1, $                ;对应Advanced Setting中的 Use Tiled Processing 1-Yes;0-No
    tile_size = 400.0, $
    wavelength_units = 'micron', $
    lambda = wavelength, $
    fwhm = fwhm, $
    input_scale = input_scale

  ;开始执行FLAASH
  flaash_obj->processImage

  call_procedure, 'envi_file_mng', id=tfid, /remove, /delete
  obj_destroy, flaash_obj
  if ~file_test(reffile_) then begin
    error='Flaash大气校正失败!'
    return,0
  endif

  call_procedure, 'envi_open_file', reffile_, r_fid=tfid
  call_procedure, 'envi_file_query', tfid, dims=dims, ns=ns, nl=nl, nb=nb
  call_procedure, 'envi_output_to_external_format', $
    fid=tfid, $
    dims=dims, $
    pos=indgen(nb), $
    out_name=reffile, $
    /tiff
  call_procedure,'envi_file_mng', id=tfid, /remove, /delete
  if ~file_test(reffile) then begin
    error='大气校正结果转TIFF格式失败!'
    return,0
  endif
  
  return,1
end

pro Rpb2rpc,rpbfile,rpcfile
  compile_opt idl2
  
  coeff = DBLARR(90)
  
  openr,lun,rpbfile,/get_lun
  line = ''
  while ~eof(lun) do begin
    readf,lun,line
    lines = strsplit(line,' =);',/extract)
    print,lines[0],'==',lines
    case strtrim(lines[0],2) of
      'lineOffset': coeff[0]=double(lines[1])
      'sampOffset': coeff[1]=double(lines[1])
      'latOffset': coeff[2]=double(lines[1])
      'longOffset': coeff[3]=double(lines[1])
      'heightOffset': coeff[4]=double(lines[1])
      'lineScale': coeff[5]=double(lines[1])
      'sampScale': coeff[6]=double(lines[1])
      'latScale': coeff[7]=double(lines[1])
      'longScale': coeff[8]=double(lines[1])
      'heightScale': coeff[9]=double(lines[1])
      'lineNumCoef': begin
        for i=0,19 do begin
          readf,lun,line
          lines = strsplit(line,' ,',/extract)
          coeff[10+i] = double(lines[0])
        endfor
      end
      'lineDenCoef': begin
        for i=0,19 do begin
          readf,lun,line
          lines = strsplit(line,' ,',/extract)
          coeff[30+i] = double(lines[0])
        endfor
      end
      'sampNumCoef': begin
        for i=0,19 do begin
          readf,lun,line
          lines = strsplit(line,' ,',/extract)
          coeff[50+i] = double(lines[0])
        endfor
      end
      'sampDenCoef': begin
        for i=0,19 do begin
          readf,lun,line
          lines = strsplit(line,' ,',/extract)
          coeff[70+i] = double(lines[0])
        endfor
      end
      else:
    endcase
  endwhile
  free_lun,lun  
  
  OPENW,rpcLun,RPCFile,/get_lun
  sign = Coeff[0] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[0],format='("LINE_OFF: '+sign+'", (d0.10), " pixels")'
  sign = Coeff[1] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[1],format='("SAMP_OFF: '+sign+'", (d0.10), " pixels")'
  sign = Coeff[2] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[2],format='("LAT_OFF: '+sign+'", (d0.10), " degrees")'
  sign = Coeff[3] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[3],format='("LONG_OFF: '+sign+'", (d0.10), " degrees")'
  sign = Coeff[4] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[4],format='("HEIGHT_OFF: '+sign+'", (d0.10), " meters")'
  sign = Coeff[5] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[5],format='("LINE_SCALE: '+sign+'", (d0.10), " pixels")'
  sign = Coeff[6] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[6],format='("SAMP_SCALE: '+sign+'", (d0.10), " pixels")'
  sign = Coeff[7] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[7],format='("LAT_SCALE: '+sign+'", (d0.10), " degrees")'
  sign = Coeff[8] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[8],format='("LONG_SCALE: '+sign+'", (d0.10), " degrees")'
  sign = Coeff[9] GE 0 ? '+' : ''
  PRINTF,rpcLun,Coeff[9],format='("HEIGHT_SCALE: '+sign+'", (d0.10), " meters")'
  
  ;LINE NUM
  FOR i=10,29 DO BEGIN
    sign = Coeff[i] GE 0 ? '+' : ''
    PRINTF,rpcLun,Coeff[i],format='("LINE_NUM_COEFF_'+STRTRIM(STRING(i+1-10),2)+': '+sign+'", E0.20)'
  ENDFOR
  ;LINE DEN
  FOR i=30,49 DO BEGIN
    sign = Coeff[i] GE 0 ? '+' : ''
    PRINTF,rpcLun,Coeff[i],format='("LINE_DEN_COEFF_'+STRTRIM(STRING(i+1-30),2)+': '+sign+'", E0.20)'
  ENDFOR
  ;SAMPLE NUM
  FOR i=50,69 DO BEGIN
    sign = Coeff[i] GE 0 ? '+' : ''
    PRINTF,rpcLun,Coeff[i],format='("SAMP_NUM_COEFF_'+STRTRIM(STRING(i+1-50),2)+': '+sign+'", E0.20)'
  ENDFOR
  ;SAMPLE DEN
  FOR i=70,89 DO BEGIN
    sign = Coeff[i] GE 0 ? '+' : ''
    PRINTF,rpcLun,Coeff[i],format='("SAMP_DEN_COEFF_'+STRTRIM(STRING(i+1-70),2)+': '+sign+'", E0.20)'
  ENDFOR
  
  FREE_LUN,rpcLun
end


function Getrpc, rpcfile
  compile_opt idl2

  line = ''
  openr,lun,rpcfile,/get_lun
  strs = strarr(90)
  readf,lun,strs
  free_lun,lun

  tempCoeff = dblarr(90)
  for i=0,89 do begin
    tempCoeff[i] = double((strsplit(strs[i],': ',/extract))[1])
  endfor

  Coeff = dblarr(93)
  Coeff[0] = tempCoeff[0]
  Coeff[1] = tempCoeff[1]
  Coeff[2] = tempCoeff[2]
  Coeff[3] = tempCoeff[3]
  Coeff[4] = tempCoeff[4]
  Coeff[5] = tempCoeff[5]
  Coeff[6] = tempCoeff[6]
  Coeff[7] = tempCoeff[7]
  Coeff[8] = tempCoeff[8]
  Coeff[9] = tempCoeff[9]
  Coeff[10:29] = tempCoeff[10:29]
  Coeff[30:49] = tempCoeff[30:49]
  Coeff[50:69] = tempCoeff[50:69]
  Coeff[70:89] = tempCoeff[70:89]
  Coeff[90:92] = [0,0,1]

  return,Coeff
end


function Orthorpc,infile,rpbfile,satellite,sensor,geofile,error=error
  compile_opt idl2
  
  if ~file_test(file_dirname(geofile)) then file_mkdir,file_dirname(geofile)

  rpcfile=file_dirname(infile)+'\'+file_basename(infile,'.tiff')+'.rpb'
  file_copy,rpbfile,rpcfile,/overwrite
  rpcfile_=file_dirname(infile)+'\'+file_basename(infile,'.tiff')+'_.rpc'
  Rpb2rpc, rpcfile, rpcfile_
  coeff=Getrpc(rpcfile_)

  zone = fix(fix(coeff[3])/6.)+31
  file_delete, rpcfile_

  case satellite of
    'GF1': begin
      if strmid(sensor,0,3) eq 'PMS' then res=8 else res=16
    end
    'GF2': res=4
    else: begin
      error='卫星类型输入有误!'
      return,0
    endelse
  endcase

  resamplemethod='bilinear';'near';'cubic';

  ; GDAL正射校正
  ;  currentdir=routine_filepath('GF_Orthorpc__define')
  ;  pricurrentdir=file_dirname(file_dirname(currentdir))
  ;  gdalwarp=pricurrentdir+'\Resource\gdalwin32-1.6\bin\'+'gdalwarp.exe'
  ;  spawn, gdalwarp + ' -of GTiff' + ' -t_srs "+proj=utm +datum=WGS84 +zone=' + strtrim(string(zone),2) + '"' + $
  ;    ' -multi' + ' -r' + ' ' + resamplemethod + ' -rpc' + ' -to' + ' "RPC_DEM=' + self.demfile + '"' + ' -tr' + $
  ;    ' ' + strtrim(string(res),2) + ' ' + strtrim(string(res),2) + ' ' + infile + ' ' + geofile, /noshell, /hide

  ; ENVI正射校正
  ;  reffile_dat=file_dirname(geofile)+'\'+file_basename(geofile,'.tiff')+'.dat'
  
  demfile=!resourcedir+'GMTED2010.jp2'
  
  envi_=call_function('envi',/headless)
  raster=envi_->OpenRaster(infile)
  dem=envi_->OpenRaster(demfile)

  tempfile=!tempdir+file_basename(geofile,'.tiff')+'_.dat'
  task=ENVITask('RPCOrthorectification')
  task.input_raster=raster
  task.dem_raster=dem
  task.output_pixel_size=res
  task.grid_spacing=10
  task.output_raster_uri=tempfile

  task->execute, error=error_

  raster->close
  dem->close
  if ~file_test(tempfile) then begin
    error=error_
    return,0
  endif

  call_procedure, 'envi_open_file', tempfile, r_fid=tfid
  call_procedure, 'envi_file_query', tfid, dims=dims, ns=ns, nl=nl, nb=nb
  call_procedure, 'envi_output_to_external_format', $
    fid=tfid, $
    dims=dims, $
    pos=indgen(nb), $
    out_name=geofile, $
    /tiff

  if ~file_test(geofile) then begin
    error='转TIFF格式失败!'
    return,0
  endif
  if file_test(rpcfile) then file_delete, rpcfile
  tfwfile=file_dirname(geofile)+path_sep()+file_basename(geofile,'.tiff')+'.tfw'
  if file_test(tfwfile) then file_delete,tfwfile
  call_procedure,'envi_file_mng',id=tfid,/remove,/delete
  call_procedure,'envi_open_file',infile,r_fid=tfid
  call_procedure,'envi_file_mng',id=tfid,/remove, /delete
  
  return,1
end