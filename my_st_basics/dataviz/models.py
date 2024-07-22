from django.db import models

# Create your models here.

class Gene(models.Model):
    barcode_id = models.CharField(max_length=200)
    global_x = models.FloatField()
    global_y = models.FloatField()
    global_z = models.FloatField()
    x = models.FloatField()
    y = models.FloatField()
    fov = models.PositiveIntegerField()
    gene = models.CharField(max_length=200)
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)
    def __str__(self):
        return f"{self.barcode_id}--{self.Gene}"
