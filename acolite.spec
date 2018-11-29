# -*- mode: python -*-

block_cipher = None

added_files = []

hiddenimports=['scipy._lib.messagestream','pandas._libs.tslibs.timedeltas']
hiddenimports+=['cftime']
hiddenimports+=['pandas._libs.tslibs.np_datetime','pandas._libs.tslibs.nattype','pandas._libs.skiplist']
hiddenimports+=['pywt._extensions._cwt']

a = Analysis(['launch_acolite.py'],
             #pathex=['/storage/Python/acolite'],
             binaries=[],
             datas=added_files,
             hiddenimports=hiddenimports,
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='acolite',
          debug=False,
          strip=False,
          upx=True,
          console=True)

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='acolite')
