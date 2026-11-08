IQMolX: IQmol with xtb calculations (see the README for what is new).

**macOS** (`IQMolX-<version>-macOS-arm64.dmg`): Apple Silicon Macs (M1 or later), macOS 15 or newer. Open the disk image and drag IQMolX to Applications. The app is not notarized by Apple: on first start right-click it and choose Open, or run `xattr -dr com.apple.quarantine /Applications/IQMolX.app`.

**Windows** (`IQMolX-<version>-windows-x64.zip`, experimental): Windows 10 or 11, 64-bit Intel/AMD. Unzip it and run `IQMolX\bin\IQMolX.exe`. Windows will warn about an unknown publisher: choose More info, then Run anyway. It compiles and runs its tests on GitHub's machines and has been tried on one real Windows 11 laptop (xtb, orbitals, frequencies, the Python engines), but it has had little testing otherwise, so please report problems.
