basefolder = [[D:\MPHYS\Data\]]
dataFileName = [[SamplePacks\]]
outFolder = [[Calcifications\]]
xdrFolder = outFolder .. [[masks\]]--string.gsub(dataFileName, [[\]], [[Xdr\]])

local function scandir(directory)
  local i, t, popen = 0, {}, io.popen
  for filename in popen('dir "'..directory..'" /o:n /b /a-d'):lines() do
    i = i + 1
    t[i] = filename
  end
  return t
end

local function median( t )
  local temp={}

  -- deep copy table so that when we sort it, the original is unchanged
  -- also weed out any non numbers
  for k,v in pairs(t) do
    if type(v) == 'number' then
      table.insert( temp, v )
    end
  end

  table.sort( temp )

  -- If we have an even number of table elements or odd.
  if math.fmod(#temp,2) == 0 then
    -- return mean value of middle two elements
    return ( temp[#temp/2] + temp[(#temp/2)+1] ) / 2
  else
    -- return middle element
    return temp[math.ceil(#temp/2)]
  end
end

local fNames = scandir(basefolder .. xdrFolder)
local medianVols = {}
local mVols = {}

for i=1,#fNames do
  print(i)
  local currentxdr = fNames[i]
  local bubs = field:new()
  --bubs:read_xdr(basefolder .. xdrFolder .. currentxdr)
  AVS:READ_XDR(bubs,basefolder .. xdrFolder .. currentxdr)
  local allVols = {}
  for j=1,bubs:max().value do
    local currentBub  = Scan:new()
    AVS:FIELD_THRESHOLD(bubs,currentBub.Data,j,j)
    allVols[j]=currentBub:volume().value
  end
  medianVols[i]=median(allVols)
  mVols[i]=medianVols[i]*bubs:max().value
end

local volFile = io.open(basefolder .. outFolder .. [[mVols.csv]], 'w')
volFile:write("ID,mVol \n")

for i=1,#fNames do
  volFile:write(fNames[i]:gsub(".xdr","") .. ',' .. mVols[i] .. '\n')
end

io.close(volFile)

print("finish")