from django.shortcuts import render
from .form import CsvModelForm
from .models import Csv
import csv
# Create your views here.

def upload_file_view(request):
    form = CsvModelForm(request.POST or None, request.FILES or None)
    if form.is_valid():
        form.save()
        form = CsvModelForm()
        print(form)
        obj = Csv.objects.filter(activated=False)[:1].get()
        print(obj)
        print("DEBUG")
        with open(obj.file_name.path, 'r') as f:
            reader = csv.reader(f)
            for i, row in enumerate(reader):
                if i==0:
                    continue
                print(row)
    return render(request,
                  'csvs/upload.html',
                  {'form': form})

