# Third-party notices and attribution

## Tectonic Explorer (Concord Consortium) — MIT

The geological model in `core/` is derived from
[Tectonic Explorer](https://github.com/concord-consortium/tectonic-explorer)
by the Concord Consortium, used under the MIT License.

> MIT License
>
> Copyright (c) 2018 Concord Consortium
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

**What was carried over:** the rock-type taxonomy and its numeric enum values
(deliberately preserved so saved scenarios remain readable across both
implementations), the crust-column and isostasy model, the subduction-vs-orogeny
decision rules, the age-as-distance-travelled convention, and the overall
phase structure of a simulation step.

**What was not:** the geodesic sphere and all its supporting machinery, the
physics integrators, the MobX/Web Worker architecture, and the Gauss–Seidel
update semantics. See `ARCHITECTURE.md` for why each was replaced.

## Unity.Mathematics — Unity Companion License

`core/` uses `Unity.Mathematics` types (`float2`, `float3`, `math`).

**This repository does not redistribute Unity.Mathematics.** It is licensed
under the [Unity Companion License](https://unity3d.com/legal/licenses/unity_companion_license),
which covers use in Unity-dependent projects but is not a permissive
redistribution licence like MIT.

- **In Unity:** add the `com.unity.mathematics` package; the `Tectonic.Core`
  assembly definition references it.
- **For the headless build:** run `scripts/fetch-deps.sh`, which clones it at a
  pinned tag into `thirdparty/` (git-ignored).
