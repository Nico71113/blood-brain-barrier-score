'use strict';
const $ = id => document.getElementById(id);
let mode = 'smiles', examples = [], current = null, batch = [], page = 0, svgUrl = null, busy = false;
const fields = ['AroR', 'HA', 'MW', 'HBA', 'HBD', 'TPSA'];
const units = {MW: 'g/mol', TPSA: 'Å²'};
const num = (v, digits=4) => typeof v === 'number' && Number.isFinite(v) ? v.toFixed(digits) : '—';
const scoreStatus = r => r.score_status === 'ok' ? 'Calculated' : r.score_status === 'invalid_input' ? 'Check input' : 'Missing input';

async function request(url, body) {
  const response = await fetch(url, {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(body)});
  if (!response.ok) { const data = await response.json(); throw new Error(data.error || 'Request failed'); }
  return response;
}
function message(text='') { $('message').textContent = text; $('message').hidden = !text; }
function loading(on) {
  busy = on;
  $('calculate').disabled = on; $('batch-calculate').disabled = on;
  document.querySelectorAll('.tab, #examples button, #local-names button').forEach(b => b.disabled = on);
  document.querySelectorAll('#single-form input, #single-form textarea, #batch-form input').forEach(input => input.disabled=on);
  $('pka').disabled=on || $('neutral').checked;
  $('calculate').textContent = on ? 'Calculating…' : 'Calculate ↗';
  $('batch-calculate').textContent = on ? 'Processing…' : 'Calculate batch ↗';
}
function invalidate() {
  if (current) { $('result-state').textContent = 'Input changed'; $('score-note').textContent = 'Recalculate to update.'; $('score-note').hidden = false; }
  $('export-single').disabled = true; $('download-svg').disabled = true; $('download-mol').disabled = true;
  $('example-note').textContent = '';
}
function switchMode(value) {
  mode = value; message();
  document.querySelectorAll('.tab').forEach(b => { const active=b.dataset.mode===mode; b.classList.toggle('active',active); b.setAttribute('aria-selected',String(active)); });
  $('single-form').hidden=mode==='batch'; $('batch-form').hidden=mode!=='batch';
  $('smiles-fields').hidden=mode!=='smiles'; $('descriptor-fields').hidden=mode!=='descriptors';
  $('local-fields').hidden=mode!=='local'; $('compound-id').required=mode==='local';
  $('name-optional').textContent=mode==='local' ? 'required' : 'optional';
  fields.forEach(key => $(key).required=mode==='descriptors');
  $('smiles').required=mode==='smiles';
  invalidate();
}
document.querySelectorAll('.tab').forEach(b => b.addEventListener('click',()=>switchMode(b.dataset.mode)));
$('single-form').addEventListener('input',invalidate);
$('neutral').addEventListener('change',()=>{ $('pka').disabled=$('neutral').checked; if($('neutral').checked) $('pka').value=''; });
function inputData() {
  const row = {compound_id:$('compound-id').value.trim(), pKa:$('pka').value || null, neutral:$('neutral').checked};
  if(mode==='smiles') row.smiles=$('smiles').value.trim();
  else if(mode==='descriptors') fields.forEach(key => row[key]=$(key).value || null);
  return row;
}
function showResult(result) {
  current=result;
  $('result-name').textContent=result.compound_id;
  $('result-state').textContent=scoreStatus(result);
  $('score-value').textContent=num(result.bbb_score);
  $('score-note').textContent=result.score_error || '';
  $('score-note').hidden=!result.score_error;
  const s=result.structure;
  if(svgUrl) { URL.revokeObjectURL(svgUrl); svgUrl=null; }
  const view=$('structure-view'); view.replaceChildren();
  if(s?.svg) {
    svgUrl=URL.createObjectURL(new Blob([s.svg],{type:'image/svg+xml'}));
    const img=document.createElement('img'); img.src=svgUrl; img.alt=`${result.compound_id} 2D structure`; view.append(img);
  } else { const p=document.createElement('p'); p.textContent=result.structure_error || 'SMILES required for 2D depiction.'; view.append(p); }
  $('formula').textContent=s?.formula ? `${s.formula} · from SMILES` : 'SMILES required';
  $('structure-source').hidden=!result.structure_source;
  $('structure-source').textContent=result.structure_source ? `Local match · ${result.structure_source} · ${result.structure_source_locator || ''}` : '';
  $('canonical-detail').hidden=!s?.canonical_smiles;
  $('canonical').textContent=s?.canonical_smiles || '';
  $('download-svg').disabled=!s?.svg; $('download-mol').disabled=!s?.molblock;
  $('descriptor-panel').hidden=false;
  $('descriptor-source').textContent=result.descriptor_source==='rdkit' ? 'RDKit' : 'Supplied';
  $('descriptor-values').replaceChildren();
  [...fields,'pKa'].forEach(key=>{
    const cell=document.createElement('div'), label=document.createElement('span'), value=document.createElement('strong');
    label.textContent=key+(units[key] ? ` (${units[key]})` : '');
    const raw=result.descriptors[key];
    value.textContent=raw==null || raw==='' ? '—' : typeof raw==='number' ? (['MW','TPSA','pKa'].includes(key) ? raw.toFixed(2) : String(raw)) : String(raw);
    cell.append(label,value); $('descriptor-values').append(cell);
  });
  $('contributions').replaceChildren();
  result.contributions.forEach(c=>{
    const row=document.createElement('div'); row.className='contribution';
    const label=document.createElement('span'); label.textContent=c.name;
    const bar=document.createElement('div'); bar.className='bar';
    const fill=document.createElement('div'); fill.className='bar-fill'; fill.style.width=`${Math.max(0,Math.min(100,c.p*100))}%`; bar.append(fill);
    const value=document.createElement('span'); value.className='number'; value.textContent=`${num(c.contribution,3)} / ${c.weight}`;
    row.append(label,bar,value); $('contributions').append(row);
  });
  $('warnings').textContent=result.warnings.join('\n');
  $('export-single').disabled=false;
}
$('single-form').addEventListener('submit',async e=>{
  e.preventDefault(); message(); loading(true);
  try { showResult(await (await request('/api/analyze',inputData())).json()); }
  catch(error) { message(error.message); }
  finally { loading(false); }
});
function download(endpoint,payload) {
  // Use an HTTP attachment; embedded WebKit browsers may not save blob URLs.
  const form=document.createElement('form'), field=document.createElement('input');
  form.method='POST'; form.action=endpoint; form.hidden=true;
  field.type='hidden'; field.name='payload'; field.value=JSON.stringify(payload);
  form.append(field); document.body.append(form); form.submit(); form.remove();
}
$('download-svg').addEventListener('click',()=>download('/api/download',{input:current.input,format:'svg'}));
$('download-mol').addEventListener('click',()=>download('/api/download',{input:current.input,format:'mol'}));
$('export-single').addEventListener('click',()=>download('/api/export',{results:[current]}));
$('export-batch').addEventListener('click',()=>download('/api/export',{results:batch}));
$('file').addEventListener('change',()=>{ $('file-name').textContent=$('file').files[0]?.name || 'No file selected'; message(); });
function readBase64(file) {
  return new Promise((resolve,reject)=>{ const reader=new FileReader(); reader.onload=()=>resolve(String(reader.result).split(',')[1]); reader.onerror=()=>reject(new Error('Cannot read file')); reader.readAsDataURL(file); });
}
$('batch-form').addEventListener('submit',async e=>{
  e.preventDefault(); const file=$('file').files[0];
  if(!file) { message('Choose a CSV or XLSX file.'); return; }
  if(file.size>16*1024*1024) { message('File exceeds 16 MB. Split it before importing.'); return; }
  message(); loading(true);
  try {
    const payload={filename:file.name,content_base64:await readBase64(file)};
    const data=await (await request('/api/batch',payload)).json(); batch=data.results; page=0; $('search').value='';
    $('structure-filter').value='all';
    $('batch-results').hidden=false;
    const s=data.summary; $('batch-summary').textContent=`${s.total} rows · ${s.scored} scored · ${s.needs_input} to review`;
    $('structure-filter').options[1].textContent=`With structure (${s.with_structure})`;
    $('structure-filter').options[2].textContent=`Without structure (${s.total-s.with_structure})`;
    renderTable();
    if(batch.length) showResult(await (await request('/api/analyze',batch[0].input)).json());
    $('batch-results').scrollIntoView({behavior:'smooth',block:'start'});
  } catch(error) { message(error.message); }
  finally { loading(false); }
});
function filteredRows() {
  const q=$('search').value.trim().toLowerCase(), state=$('structure-filter').value;
  return batch.filter(r=>r.compound_id.toLowerCase().includes(q) && (state==='all' || (r.structure_status==='ok')===(state==='available')));
}
function renderTable() {
  const filtered=filteredRows(), pages=Math.max(1,Math.ceil(filtered.length/50)); page=Math.max(0,Math.min(page,pages-1));
  $('batch-body').replaceChildren();
  filtered.slice(page*50,page*50+50).forEach(r=>{
    const tr=document.createElement('tr'); tr.tabIndex=0;
    const structureLabel=r.structure_status==='ok' ? (r.structure_origin==='local_association' ? 'Matched' : 'Available') : (r.structure_status==='invalid_structure' ? 'Invalid SMILES' : 'Missing SMILES');
    const values=[r.row_id,r.compound_id,num(r.bbb_score),r.bbb_score_input ?? '—',structureLabel,r.score_error || r.structure_error || (r.warnings.length ? r.warnings.join('; ') : 'OK')];
    values.forEach((v,i)=>{const td=document.createElement('td'); td.textContent=String(v); if(i===2||i===3)td.className='numeric'; if(i===5)td.className='note'; tr.append(td);});
    const select=async()=>{
      if(busy) return;
      loading(true); message();
      try { const detail=await (await request('/api/analyze',r.input)).json(); detail.row_id=r.row_id; showResult(detail); $('result-name').scrollIntoView({behavior:'smooth',block:'center'}); }
      catch(error){message(error.message);} finally{loading(false);}
    };
    tr.addEventListener('click',select); tr.addEventListener('keydown',e=>{if(e.key==='Enter'){e.preventDefault();select();}}); $('batch-body').append(tr);
  });
  $('page-info').textContent=`Page ${page+1} / ${pages} · ${filtered.length} ${filtered.length===1 ? 'row' : 'rows'}`;
  $('prev-page').disabled=page===0; $('next-page').disabled=page>=pages-1;
}
$('search').addEventListener('input',()=>{page=0;renderTable();});
$('structure-filter').addEventListener('change',()=>{page=0;renderTable();});
$('prev-page').addEventListener('click',()=>{page--;renderTable();});
$('next-page').addEventListener('click',()=>{page++;renderTable();});
async function loadExample(i) {
  if(busy) return;
  const e=examples[i]; switchMode('smiles');
  $('compound-id').value=e.compound_id; $('smiles').value=e.smiles; $('pka').value=e.pKa ?? ''; $('neutral').checked=false; $('pka').disabled=false;
  $('example-note').textContent=e.note; message(); loading(true);
  try { showResult(await (await request('/api/analyze',e)).json()); }
  catch(error) { message(error.message); } finally { loading(false); }
}
async function init() {
  try {
    const indexResponse=await fetch('/api/structures-status');
    if(!indexResponse.ok) throw new Error('Cannot read the local structure index.');
    const index=await indexResponse.json();
    document.querySelector('[data-mode="local"]').hidden=index.records===0;
    $('library-status').hidden=!index.records;
    $('library-status').textContent=`${index.records} verified local structures available.`;
    index.names.forEach(name=>{
      const b=document.createElement('button');b.type='button';b.textContent=name;
      b.addEventListener('click',async()=>{
        if(busy) return;
        $('compound-id').value=name; $('pka').value=''; $('neutral').checked=false; $('pka').disabled=false;
        loading(true);message();
        try {showResult(await (await request('/api/analyze',{compound_id:name})).json());}
        catch(error){message(error.message);}finally{loading(false);}
      });
      $('local-names').append(b);
    });
    const r=await fetch('/api/examples'); if(!r.ok) throw new Error('Cannot load examples.'); examples=await r.json();
    examples.forEach((e,i)=>{const b=document.createElement('button');b.type='button';b.textContent=e.compound_id;b.addEventListener('click',()=>loadExample(i));$('examples').append(b);});
  }catch(error){message(error.message);}
}
init();
